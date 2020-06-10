#### Time Scale Analysis ################################################################

"""
Calculate threshold
Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation 
"""
calculate_threshold(vm_arr::AbstractArray where T) = sum(vm_arr)/length(vm_arr) + 4*std(vm_arr)

"""
This function acts to calculate the distance between points in a single BitArray. 
Very rarely is the first point part of a spike (in which case there is a fallback), 
and because of this clip is set to remove the first interval. 
"""
function count_intervals(spike_trace::BitArray{1}; clip = 2)
    isi = Array{Float64}([])
    count = 0
    
    #In the off chance that the first spike occurs right as the first time point, then clipping is cancelled. 
    if spike_trace[1] == 1
        clip = 1
    end
    
    for spike in spike_trace
        if spike == 0
            count += 1
        elseif spike == 1 && count == 0
            count = 0
        elseif spike == 1 && count > 0
            push!(isi, count)
            count = 0
        end
    end
    isi[clip:end]
end

function count_intervals(spike_arr::BitArray{2}, clip = 2)
    isi_arr = Array{Float64}([])
    count = 0
    for idx_y = 1:size(spike_arr,1)
        all_intervals = extract_interval(spike_arr[idx_y, :], clip = 2)
        push!(isi_arr, all_intervals...)
    end
    isi_arr
end

function count_intervals(spike_arr::BitArray{3})
    isi_arr = Array{Float64}([])
    count = 0
    for idx_y = 1:size(spike_arr,1)
        for idx_x = 1:size(spike_arr, 2)
            all_intervals = extract_interval(spike_arr[idx_y, idx_x, :], clip = 2)
            push!(isi_arr, all_intervals...)
        end
    end
    isi_arr
end

"""
This function returns all the time stamps in a spike or burst array
"""
function get_timestamps(spike_array::BitArray{1}; dt = 1.0)
    intervals = count_intervals(spike_array) .* dt
    durations = count_intervals(spike_array .!= 1.0) .* dt
    first_point = findfirst(x -> x==1, spike_array)-1|>Float64
    current_point = first_point
    points = Tuple[(current_point, current_point+durations[1])]
    
    for idx = 1:length(intervals)
        current_point += intervals[idx] + durations[idx+1]
        push!(points, (current_point-durations[idx+1], current_point))
    end
    return points
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
"""
function max_interval_algorithim(spike_array::BitArray{1}; ISIstart = 500, ISIend = 500, IBImin = 1000, DURmin = 100, SPBmin = 4, dt = 1.0, verbose = false)
    burst_timestamps = Array{Tuple,1}([])
    DUR_list = Array{Float64,1}([])
    SPB_list = Array{Float64,1}([])
    IBI_list = Array{Float64,1}([])
    
    intervals = count_intervals(spike_array) .* dt
    timestamps = get_timestamps(spike_array; dt = dt)
    bursting = false
    burst_start = 0.0
    burst_end = 0.0
    IBI = 0.0
    SPB = 0
    idx = 1
    for i = 1:length(intervals)
        if bursting == false && intervals[i] <= ISIstart
            bursting = true
            burst_start = timestamps[i][1]
        elseif bursting == true && intervals[i] >= ISIend
            bursting = false
            burst_end = timestamps[i][2]
            IBI = intervals[i]
            DUR = (burst_end - burst_start)
            if IBI >= IBImin && DUR >= DURmin && SPB >= SPBmin
                if verbose
                    println("Timestamp $idx: $burst_start -> $burst_end, DUR $idx: $DUR,  SPB $idx: $SPB, IBI $idx: $IBI,")
                end
                push!(burst_timestamps, (burst_start, burst_end))
                push!(DUR_list, DUR)
                push!(SPB_list, SPB)
                push!(IBI_list, IBI)
                SPB = 0
                idx+=1
            end  
        end
        if bursting == true
            SPB += 1
        end
    end
    #This algorithim usually leaves one last burst off because it has no end point. We can add this
    DUR = (timestamps[end][2] - burst_start)
    if DUR >= DURmin && SPB >= SPBmin && bursting = true
        if verbose
            println("Timestamp  $idx: $burst_start -> $(timestamps[end][2]), DUR $idx: $DUR, SPB $idx: $SPB, IBI $idx: Unknown")
        end
        push!(burst_timestamps, (burst_start, timestamps[end][2]))
        push!(DUR_list, DUR)
        push!(SPB_list, SPB)
    end
    burst_timestamps, DUR_list, SPB_list, IBI_list
end

function timescale_analysis(vm_trace; dt = 10.0, verbose = false)
    sim_thresh = calculate_threshold(vm_trace)
    spike_array = (vm_trace .> sim_thresh);
    if !any(spike_array .== 1.0)
        if verbose
            println("No spikes detected")
        end
        return fill(NaN, 6)
    else
        timestamps = get_timestamps(spike_array; dt = dt);
        durations = map(x -> x[2]-x[1], timestamps)
        avg_spike_dur = sum(durations)/length(durations)
        std_spike_dur = std(durations)

        burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array; verbose = verbose, dt = dt);
        avg_burst_dur = sum(dur_list)/length(dur_list);
        std_burst_dur = std(dur_list);

        avg_ibi = sum(ibi_list)/length(ibi_list)
        std_ibi = std(ibi_list)
        if verbose
            println("The average spike duration is $(round(avg_spike_dur, digits = 2)) ms +- $(round(std_spike_dur, digits = 2)) ms")
            println("The average burst duration is $(round(avg_burst_dur;digits = 2)) ms +- $(round(std_burst_dur;digits = 2)) ms")
            println("The average interburst interval is $(round(avg_ibi;digits = 2)) ms +- $(round(std_ibi;digits = 2)) ms")
        end
        return [avg_spike_dur, std_spike_dur, avg_burst_dur, std_burst_dur, avg_ibi, std_ibi]
    end
end


function calculate_STTC(signal1::BitArray{1}, signal2::BitArray{1}; Δt::Float64 = 50.0, dt = 1.0)
    n_spike1 = sum(signal1);
    n_spike2 = sum(signal2);
    burst_conv1, burst_idxs1 = convolve_bursts(signal1; θr = Δt, dt = dt, include_theta = true);
    burst_conv2, burst_idxs2 = convolve_bursts(signal2; θr = Δt, dt = dt, include_theta = true);
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

function calculate_STTC(spike_array::BitArray{2}; Δt::Float64 = 50.0, dt::Float64 = 1.0)
    n_traces = size(spike_array, 1)
    corr_matrix = zeros(n_traces, n_traces)
    for i = 1:n_traces
        for j = 1:n_traces
            corr_matrix[i,j] = calculate_STTC(spike_array[i, :], spike_array[j, :])
        end
    end
    corr_matrix
end








############# Extracting Wave fronts #######################################################
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
    burst_arr = convolve_bursts(spike_arr)
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
This function will be replaced with max_interval_algorithim, but for now is a placeholder for all other functions still utilizing it. 
"""
function convolve_bursts(spike_arr::BitArray{1}; θr = 500.0, IBI = 1000.0, min_dur = 100.0, dt = 10.0, include_theta = false)
    spike_timer = 0.0
    duration_timer = 0.0
    burst_start = 0
    burst_end = 0
    burst_inds = Array{Tuple}([])
    burst_arr = zeros(Int, length(spike_arr))
    for idx in 1:length(spike_arr)
        if spike_arr[idx] == 1 && duration_timer == 0
            spike_timer = θr
            if burst_start == 0
                burst_start = idx
            else
                burst_end = idx
            end
        elseif spike_arr[idx] == 0 && spike_timer > 0
            spike_timer -= dt
        else
            spike_timer = 0.0
            
            #If a burst has concluded and duration is greater than 0, countdown the duration timer
            if duration_timer > 0
                duration_timer -= dt
            end
            
            if burst_start > 0 && burst_end != 0 && burst_end - burst_start > min_dur
                
                push!(burst_inds, (burst_start, burst_end))
                if include_theta == true
                    burst_start = max(1, round(Int, burst_start-(θr/dt)))
                    burst_end = round(Int, min(length(burst_arr), burst_end+(θr/dt)))
                end                    
                burst_arr[burst_start:burst_end] .= 1
                #This signifies that a burst has ended and we can begin a duration timer
                duration_timer = IBI
            end
            burst_start = 0
            burst_end = 0
        end
    end
    burst_arr, burst_inds
end

function convolve_bursts(spike_arr::BitArray{2}; θr = 500.0,  IBI = 1000.0, min_dur = 100.0, dt = 10.0, include_theta = false)
    ret_arr = similar(spike_arr)
    burst_inds = Array{Tuple}([])
    for i = 1:size(spike_arr, 1)
        ret_arr[i,:], b_idxs = convolve_bursts(spike_arr[i, :]; θr = θr, IBI = IBI, min_dur = min_dur, dt = dt, include_theta = include_theta)
        push!(burst_inds, (i, b_idxs...))
    end
    ret_arr, burst_inds
end

function convolve_bursts(spike_arr::BitArray{3}; θr =  500.0, IBI = 1000.0, min_dur = 100.0, dt = 10.0, include_theta = false)
    ret_arr = similar(spike_arr)
    burst_inds = Array{Tuple}([])
    for i = 1:size(spike_arr, 1)
        for j = 1:size(spike_arr, 2)
            ret_arr[i,j,:], b_idxs = convolve_bursts(spike_arr[i,j, :];  θr = θr, IBI = IBI, min_dur = min_dur, dt = dt, include_theta = include_theta)
            push!(burst_inds, ((i, j), b_idxs...))
        end
    end
    ret_arr, burst_inds
end
