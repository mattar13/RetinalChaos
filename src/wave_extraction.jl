"""
Function to extract 3D points
"""
function extract_points(points::Array{CartesianIndex{3},1})
    xs = map(x -> x[1], points)
    ys = map(x -> x[2], points)
    zs = map(x -> x[3], points)
    return xs, ys, zs
end

"""
Function to extract 2D points
"""
function extract_points(points::Array{CartesianIndex{2},1})
    xs = map(x -> x[1], points)
    ys = map(x -> x[2], points)
    return xs, ys
end

"""
Calculate threshold
"""
calculate_threshold(vm_arr::AbstractArray where T) = sum(vm_arr)/length(vm_arr) + 4*std(vm_arr)


"""
This function extracts the interspike intervals in a single wave.
    It is designed to clip any intervals at the beginning (as they may or may not be part of a previous burst)
    and clip events at the end (because we don't know when the next burst will occur)
    - If interval_or_count is set to 1, then the isi extracts a spike count
"""
function extract_isi(spike_trace::BitArray{1}; dt = 1.0)
    isi = Array{Float64}([])
    count = 0
    for spike in spike_trace
        if spike == 0
            count += 1
        elseif spike == 1 && count == 0
            count = 0
        elseif spike == 1 && count > 0
            push!(isi, count*dt)
            count = 0
        end
    end
    isi
end

"""
This function calculates the interspike interval
The input must always be a thresholded spike array
"""
function calculate_isi(spike_arr::BitArray{3};  dt = 10.0)
    isi_arr = Array{Float64}([])
    count = 0
    for idx_y = 1:size(spike_arr,1)
        for idx_x = 1:size(spike_arr, 2)
            all_intervals = extract_isi(spike_arr[idx_y, idx_x, :]; dt = dt)
            push!(isi_arr, all_intervals...)
        end
    end
    sort(isi_arr)
end

"""
Function to find the maximum interspike interval
"""
function find_maxISI(spike_arr::BitArray{3}; dt = 10.0, rank = 0.75, ci = 0.95, spike_dist = Exponential)
    isi = calculate_isi(spike_arr; dt = dt) #Calculate all ISIs
    spike_intervals, _ = rank_isi(isi; dt = dt, rank = rank)#all intervals below rank 75 are spikes
    spike_fit = fit(spike_dist, spike_intervals) #Fit the spikes to a normal distribution
    quantile(spike_fit, ci) #ANything beyond the 95 percentile is not a interspike distance
end

function rank_isi(isi::Vector{Float64}; dt = 10.0, rank = 0.75)
    unique_isi = unique(isi) #Calculate all ISI ranks starting from shortests
    ranking = unique_isi[round(Int, length(unique_isi)*(1-rank))]
    spike_intervals = isi[findall(x -> x < ranking, isi)] #all intervals below rank 75 are spikes
    burst_intervals = isi[findall(x -> x >= ranking, isi)]
    return spike_intervals, burst_intervals
end

"""
This uses the minimum spike calculated above to convolve the spike trains into bursts
"""
function convolve_bursts(spike_arr::BitArray{1}, θr; dt = 10.0, include_theta = false)
    burst_timer = 0.0
    burst_start = 0
    burst_end = 0
    burst_inds = Array{Tuple}([])
    burst_arr = zeros(Int, length(spike_arr))
    for idx in 1:length(spike_arr)
        if spike_arr[idx] == 1
            burst_timer = θr
            if burst_start == 0
                burst_start = idx
            else
                burst_end = idx
            end
        elseif spike_arr[idx] == 0 && burst_timer > 0
            burst_timer -= dt
        else
            burst_timer = 0.0
            if burst_start > 0 && burst_end != 0
                push!(burst_inds, (burst_start, burst_end))
                if include_theta == true
                    burst_start = max(0, round(Int, burst_start-(θr/dt)))
                    burst_end = round(Int, min(length(burst_arr), burst_end+(θr/dt)))
                end                    
                burst_arr[burst_start:burst_end] .= 1
            elseif burst_start > 0 && burst_end == 0 
                #push!(burst_inds, (burst_start, burst_start))
                #burst_arr[burst_start] = 1
            end
            burst_start = 0
            burst_end = 0
        end
    end
    burst_arr, burst_inds
end

function convolve_bursts(spike_arr::BitArray{2}, θr; dt = 10.0, include_theta = false)
    ret_arr = similar(spike_arr)
    burst_inds = Array{Tuple}([])
    for i = 1:size(spike_arr, 1)
        ret_arr[i,:], b_idxs = convolve_bursts(spike_arr[i, :], θr; dt = dt, include_theta = include_theta)
        push!(burst_inds, (i, b_idxs...))
    end
    ret_arr, burst_inds
end

function convolve_bursts(spike_arr::BitArray{3}, θr; dt = 10.0, include_theta = false)
    ret_arr = similar(spike_arr)
    burst_inds = Array{Tuple}([])
    for i = 1:size(spike_arr, 1)
        for j = 1:size(spike_arr, 2)
            ret_arr[i,j,:], b_idxs = convolve_bursts(spike_arr[i,j, :], θr, include_theta = include_theta)
            push!(burst_inds, ((i, j), b_idxs...))
        end
    end
    ret_arr, burst_inds
end

"""
This function seperates bursts from waves It completes tasks in this order
"""
function wave_finder(burst_arr::BitArray{3} where T)
    markers = label_components(burst_arr)
    n_waves = maximum(markers)
    burst_events = Vector{Vector{CartesianIndex{3}}}([])
    wave_events = Vector{Vector{CartesianIndex{3}}}([])
    @showprogress 0.5 "Finding waves..." for wave_idx = 1:n_waves
        pᵢ = findall(x -> x == wave_idx, markers)
        (xs, ys, zs) = pᵢ |> extract_points

        if unique(xs) |> length > 1 || unique(ys) |> length > 1
            push!(wave_events, pᵢ)
        else
            push!(burst_events, pᵢ)
        end
    end
    burst_events, wave_events
end


"""
This function calculates the percentage of spikes in a burst for all Spike and Burst traces
"""
function spike_in_burst(burst_arr::BitArray{3}, spike_arr::BitArray{3})
    pib = Array{Float64}([])
    for idx_y = 1:size(burst_arr, 1)
        for idx_x = 1:size(burst_arr, 2)
            push!(pib, spike_in_burst(burst_arr[idx_y, idx_x, :], spike_arr[idx_y, idx_x, :]))
        end
    end
    pib
end

function spike_in_burst(burst_arr::BitArray{2}, spike_arr::BitArray{2})
    pib = Array{Float64}([])
    for idx_y = 1:size(burst_arr, 1)
        push!(pib, spike_in_burst(burst_arr[idx_y, :], spike_arr[idx_y, :]))
    end
    pib
end

function spike_in_burst(burst_arr::BitArray{1}, spike_arr::BitArray{1})
    trace_in_burst = burst_arr .* spike_arr
    trace_not_burst = (burst_arr.==0.0) .* spike_arr
    1 - sum(trace_not_burst)/sum(trace_in_burst)
end


"""
These functions calculate the ISIs of all of the bursts
"""
function calculate_burstISI(spike_arr::BitArray{1}, θr; dt = 10.0)
    bs, b_ids = convolve_bursts(spike_arr, θr; dt = dt)
    burst_ISIs = Array{Float64}([])
    for (b_start, b_end) in b_ids
        if b_start < b_end
            burst_duration = (b_end * 0.01) - (b_start * 0.01)
            burst_i = spike_arr[b_start:b_end]
            n_spikes = sum(burst_i)
            burst_isi_list = extract_isi(burst_i; dt = 0.01)
            push!(burst_ISIs, burst_isi_list...)
        end
    end
    burst_ISIs
end

function calculate_burstISI(spike_arr::BitArray{2}, θr; dt = 10.0)
    burst_ISIs = Array{Float64}([])
    for idx_y = 1:size(spike_arr, 1)
        push!(burst_ISIs, calculate_burstISI(spike_arr[idx_y, :], θr; dt = dt)...)
    end
    burst_ISIs
end

function calculate_burstISI(spike_arr::BitArray{3}, θr; dt = 10.0)
    burst_ISIs = Array{Float64}([])
    for idx_y = 1:size(spike_arr, 1)
        for idx_x = 1:size(spike_arr, 2)
            push!(burst_ISIs, calculate_burstISI(spike_arr[idx_y, idx_x, :], θr; dt = dt)...)
        end
    end
    burst_ISIs
end


"""
This records stats according to table 2 in Xu Et al. 2016

Stats included in the function
Each stat has a mean and a standard deviation
1) Firing frequency
2) Interburst Interval
3) Burst frequency
4) Burst duration
5) % of spike in burst
6) Interspike interval during burst
7) Wave frequency
8) Wave duration
9) Wave size (% of cells)
10) Wave Interspike Interval

Other stats we could add
a) firing rate in wave
b) firing rate in bursts
"""
function run_wavestats(timestamp, vm_arr; dt = 0.01, default_θ = :none)

    threshold = calculate_threshold(vm_arr)
    spike_arr = vm_arr .>= threshold
    allISI = calculate_isi(spike_arr; dt = dt)
    #calculate the firing frequency and standard deviation
    n_FFq = length(allISI)
    FiringFrequency = sum(allISI)/n_FFq
    FiringFrequencyStd = std(allISI)
    #Calculate the interburst interval
    ISI, IBI = rank_isi(allISI)
    n_IBI = length(IBI)
    InterburstInterval = sum(IBI)/n_IBI
    InterburstIntervalStd = std(IBI)
    #Calculate the burst frequency
    n_BFq = n_IBI
    BurstFrequency = 60/InterburstInterval#in bursts per minute
    BurstFrequencyStd = 1/InterburstIntervalStd
    #Find the burst/spike threshold
    if default_θ == :none
        θr = find_maxISI(spike_arr)
    else
        θr = default_θ
    end
    burst_arr = convolve_bursts(spike_arr, θr)
    burst_length = calculate_isi(burst_arr.==0; dt = dt) #All burst lengths
    #Calculate the burst durations
    n_BD = length(burst_length)
    BurstDuration = sum(burst_length)/n_BD #Average burst duration
    BurstDurationStd = std(burst_length) #Std
    #Calculate the percentage of spikes in a burst
    spike_ib = spike_in_burst(burst_arr, spike_arr)
    n_PB = length(spike_ib)
    PercentBurstSpikes = sum(spike_ib)/n_PB
    PercentBurstSpikesStd = std(spike_ib)
    #Calculate the interspike interval within a burst
    burst_isis = calculate_burstISI(spike_arr, θr)
    n_ISI = length(burst_isis)
    InterspikeInterval = sum(burst_isis)/n_ISI
    InterspikeIntervalStd = std(burst_isis)
    #Go through the burst array and find wave events
    burst_events, wave_events = wave_finder(burst_arr)
    #Calculate the percentage of cells each wave includes
    if length(wave_events) > 0
        wave_per_timepoint = zeros(size(burst_arr,3))
        percent_cells = []; wave_durations = []; wave_isis = []
        for wave in wave_events
            (xs, ys, zs) = wave |> extract_points
            time_points = unique(zs)
            wave_per_timepoint[time_points] .= 1.0
            wave_dur = (time_points[end]*dt) - (time_points[1]*dt)
            push!(wave_durations, wave_dur)
            push!(percent_cells, length(unique(xs)) * length(unique(ys))/(size(spike_arr, 1) * size(spike_arr, 2)))
            all_wave_isis = calculate_isi(spike_arr[:,:,unique(zs)]; dt = dt)
            push!(wave_isis, all_wave_isis...)
        end
        #Calculate wave frequencies
        IWI = extract_isi(wave_per_timepoint .> 0.0; dt = dt)
        n_WFq = length(IWI)
        WaveFrequency = sum(IWI)/n_WFq
        WaveFrequencyStd = std(IWI)
        #Calculate wave durations
        n_WD = length(wave_durations)
        WaveDuration = sum(wave_durations)/n_WD
        WaveDurationStd = std(wave_durations)
        #Calculate Wave size (% of cells)
        n_WS = length(percent_cells)
        WaveSize = sum(percent_cells)/n_WS
        WaveSizeStd = std(percent_cells)
        #calculate mean wave interspike interval
        n_ISIw = length(wave_isis)
        WaveInterspikeInterval = sum(wave_isis)/n_ISIw
        WaveInterspikeIntervalStd = std(wave_isis)
    else
        println("No detectable waves")
        n_WFq = WaveFrequency = WaveFrequencyStd = 0.0
        #Calculate wave durations
        n_WD = WaveDuration = WaveDurationStd = 0.0
        #Calculate Wave size (% of cells)
        n_WS = WaveSize = WaveSizeStd = 0.0
        #calculate mean wave interspike interval
        n_ISIw = WaveInterspikeInterval = WaveInterspikeIntervalStd = 0.0
    end

    #Make a blank datasheet
    df_stats = DataFrame(
        Date = DateTime[],
        n_FFq = Int64[], FiringFrequency = Float64[], FiringFrequencyStd = Float64[],
        n_IBI = Int64[], InterburstInterval = Float64[], InterburstIntervalStd = Float64[],
        n_BFq = Int64[], BurstFrequency = Float64[], BurstFrequencyStd = Float64[],
        n_BD  = Int64[], BurstDuration = Float64[], BurstDurationStd = Float64[],
        n_PB  = Int64[], PercentBurstSpikes = Float64[], PercentBurstSpikesStd = Float64[],
        n_ISI = Int64[], InterspikeInterval = Float64[], InterspikeIntervalStd = Float64[],
        n_WFq = Int64[], WaveFrequency = Float64[], WaveFrequencyStd = Float64[],
        n_WD  = Int64[], WaveDuration = Float64[], WaveDurationStd = Float64[],
        n_WS  = Int64[], WaveSize = Float64[], WaveSizeStd = Float64[],
        n_ISIw = Int64[], WaveInterspikeInterval = Float64[], WaveInterspikeIntervalStd = Float64[]
    )
    data_vector = [
        timestamp,
        n_FFq, FiringFrequency, FiringFrequencyStd,
        n_IBI, InterburstInterval, InterburstIntervalStd,
        n_BFq , BurstFrequency, BurstFrequencyStd,
        n_BD  , BurstDuration, BurstDurationStd,
        n_PB  , PercentBurstSpikes, PercentBurstSpikesStd,
        n_ISI , InterspikeInterval, InterspikeIntervalStd,

        n_WFq , WaveFrequency, WaveFrequencyStd,
        n_WD  , WaveDuration, WaveDurationStd,
        n_WS  , WaveSize, WaveSizeStd,
        n_ISIw , WaveInterspikeInterval, WaveInterspikeIntervalStd
    ]
    push!(df_stats, data_vector)
    df_stats
end

"""
This function finds local maxima in a graph. Usually requires noise to be cleaned with a loess function
"""
function findlocalmaxima(array)
    maxima = Int64[]
    for idx = 2:length(array)-1
        if array[idx-1] < array[idx] > array[idx+1]
            push!(maxima, idx)
        end
    end
    maxima
end

"""
This method uses the Loess filter function to find peaks within a histogram or waveform
"""
function find_baselines(array)
    thresh = calculate_threshold(vm[reduced]);
    bins = collect(LinRange(minimum(vm[reduced])-1,thresh, 1000))
    h_fit = fit(Histogram, vm, bins, closed = :left);
    h_edges = Float64.(h_fit.edges[1][2:end]);
    h_weights = Float64.(h_fit.weights);
    h_model = loess(h_edges, h_weights; span = 0.15, degree = 2);
    smooth_weights = Loess.predict(h_model, h_edges);
    findlocalmaxima(smooth_weights)
end



function calculate_STTC(signal1::BitArray{1}, signal2::BitArray{1}; θr::Float64 = 50.0, dt = 1.0)
    n_spike1 = sum(signal1);
    n_spike2 = sum(signal2);
    burst_conv1, burst_idxs1 = convolve_bursts(signal1, θr; dt = dt, include_theta = true);
    burst_conv2, burst_idxs2 = convolve_bursts(signal2, θr; dt = dt, include_theta = true);
    T1 = sum(burst_conv1)/length(signal1)*dt;
    T2 = sum(burst_conv2)/length(signal2)*dt;

    s_b1 = 0.0
    for (sta, en) in burst_idxs2
        s_b1 += sum(signal1[sta:en])
    end
    P1 = s_b1/n_spike1

    s_b2 = 0.0
    for (sta, en) in burst_idxs1
        s_b2 += sum(signal2[sta:en])
    end
    P2 = s_b2/n_spike2
    1/2*((P1-T2)/(1-P1*T2) + (P2-T1)/(1-P2*T1))
end

"""
max_val = maximum(isi); min_val = minimum(isi); dBin = 10.0;
bins = collect(min_val:dBin:max_val)
isi_hist = fit(Histogram, isi, bins, closed = :left)
plot(isi_hist, xlabel = "ISI (s)", ylabel = "LogFrequency")

e_isi = isi_hist.edges[1][2:end]
w_isi = log.(isi_hist.weights)
plot(e_isi, w_isi)
model = loess(e_isi, w_isi, span = 0.33, degree = 4)
isi_loess = Loess.predict(model, e_isi)
plot!(e_isi, isi_loess*10)
pi_k = findlocalmaxima(isi_loess)

under_mcv = findall(x -> x<mcv, e_isi)
ibp_k = argmax(isi_hist.weights[under_mcv])

#now we need to find the C_peaks (will occur after the MCV)
e_isi_mcv = isi_hist.edges[1][under_mcv[end]+1:end]
w_isi_mcv = isi_hist.weights[under_mcv[end]:end]
model = loess(e_isi_mcv, w_isi_mcv, span = 0.33, degree = 4)
isi_loess = Loess.predict(model, e_isi_mcv);
pi_k = findlocalmaxima(isi_loess)


argsmin(series) = findall(x -> x == minimum(series), series)
min_idxs = argsmin(isi_hist.weights[ibp_k:pi_k[1]])
map(isi_min -> void(isi_min, c_pi[1], c_ibp), isi_hist.edges[1][min_idxs])

maxISI = isi_hist.edges[1][min_idx]

c_pi = isi_hist.weights[pi_k][1]
c_ibp = isi_hist.weights[ibp_k]
void(maxISI, c_pi[1], c_ibp)

plot(e_isi_mcv, w_isi_mcv)
vline!(e_isi_mcv[pi_k], lw = 2.0)


"""
not_good() = println("This stuff is not good")
