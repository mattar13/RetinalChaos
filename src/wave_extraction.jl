#### Time Scale Analysis ################################################################

"""
    calculate_threshold(sol::DiffEqBase.AbstractODESolution, Z::Int64)

Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation 
"""
calculate_threshold(vm_arr::AbstractArray; Z::Int64 = 4) = [sum(vm_arr)/length(vm_arr) + Z*std(vm_arr)]

function calculate_threshold(sol::DiffEqBase.AbstractODESolution{T, N, S}, rng::Tuple{Float64, Float64}; 
        idx::Int64 = 1, Z::Int64 = 4, dt::Float64 = 100.0,
    ) where {T, N, S}
    # We need to convert the dt into the correct form
    dt = convert(T, dt)
    #We want to check how many dimensions the simulation is 
    if length(size(sol.prob.u0)) == 1
        data_section = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
        #println(data_section |> size)
        return calculate_threshold(data_section; Z = Z)
    else  
        n = length(collect(rng[1]:dt:rng[2]))
        mean = sum(sol(rng[1]:dt:rng[2]), dims = 2)./n
        dev = Z * std(sol(rng[1]:dt:rng[2]))
        return (mean .+ dev) |> Array
    end
end

function calculate_threshold(sol::DiffEqBase.AbstractODESolution{T, N, S}; kwargs...) where {T, N, S} 
    t_rng = (sol.t[1]|> Float64, sol.t[end] |> Float64)
    ans = calculate_threshold(sol, t_rng; kwargs...)
    return ans
end

"""
This function returns all the time stamps in a spike or burst array
"""

function get_timestamps(spike_array::Vector{Bool}, tseries::Vector{T}) where T <: Real
    diff_vals = map(i -> (spike_array[i]-spike_array[i+1]), 1:length(spike_array)-1)
    diff_starts = findall(x -> x==-1, diff_vals).+1 #This is a list of all the starting points in the array
    diff_ends = findall(x -> x==1, diff_vals) #This is a list of all the ending points in the array
    #If there is one more end than start than we have to remove the last point from the diff_starts

    if spike_array[1] #This means we start out in a spike and will most likely end o
        #println("We started out in a spike, the first value will be an end spike")
        diff_ends = diff_ends[2:end]
    end

    if length(diff_starts) > length(diff_ends)  #This happens because an end point was cutoff
        diff_starts = diff_starts[1:length(diff_ends)]
    elseif length(diff_starts) < length(diff_ends) #This happens because a start point was cutoff
        diff_ends = diff_ends[2:end]
    end

    return hcat(tseries[diff_starts], tseries[diff_ends])
end

function get_timestamps(spike_array::Matrix{Bool}, timestamps::Vector{T}) where T <: Real
    tstamps = Vector{Matrix{T}}(undef, size(spike_array,1))
    for i in 1:size(spike_array, 1)
        tstamps[i] = get_timestamps(spike_array[i, :], timestamps)
    end
    return tstamps
end

# For if no range has been provided but a threshold has
function get_timestamps(sol::DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}; 
        idx::Int64 = 1, dt::Float64 = -Inf
    ) where T <: Real
    if dt == -Inf #Adaptive timestamp finding
        #we need to find a way to do 2D arrays vs 1D
        if length(size(sol.u)) == 1
            data_select = sol |> Array |> Array{T}
            spike_array = (data_select .> threshold) |> Array
            tstamps = get_timestamps(spike_array, sol.t |> Array)
            return tstamps
        elseif length(size(sol.u)) == 2
            data_select = sol |> Array |> Array{T}
            spike_array = (data_select .> threshold) |> Array
            tstamps = get_timestamps(spike_array, sol.t |> Array)
            return tstamps
        end
    else   
        get_timestamps(sol, threshold, (sol.t[1], sol.t[end]); idx = idx, dt = dt)
    end
end

"""
If dt is set to Inf, the algorithim acts adaptive
"""
function get_timestamps(sol::DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = Inf
    ) where T <: Real
    # Lets do this the integrated way
    # First we need to extract the spike array
    if length(size(sol.prob.u0)) == 1
        data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
        spike_array = (data_select .> threshold[1])
        get_timestamps(spike_array; dt = dt)
    else
        #For each data trace we have a corresponding threshold
        timestamps = Tuple{T, T}[]
        @showprogress "Extracting timestamps: " for c_idx in 1:size(sol,1)
            data_select = sol(rng[1]:dt:rng[2], idxs = c_idx) |> Array
            spike_array = (data_select .> threshold[c_idx])
            stamps = get_timestamps(spike_array; dt = dt)
            push!(timestamps, stamps...)
        end
        return timestamps
    end   
end

# get timestamps will not work adaptively with a region defined
function get_timestamps(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T,T}; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real
    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    get_timestamps(sol, threshold, rng; dt = dt)
end

# For if no range has been provided
function get_timestamps(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real
    threshold = calculate_threshold(sol; idx = idx, Z = Z)
    get_timestamps(sol, threshold; idx = idx, dt = dt)
end

function extract_interval(timestamps::Matrix{T}; 
        max_duration = 10e5, max_interval = 10e5, 
        min_interval = 0.0, min_duration = 0.0
    ) where T <: Real
    durations = timestamps[:, 2] .- timestamps[:,1]
    lagged_starts = timestamps[2:end,1]
    lagged_ends = timestamps[1:end-1,2]
    intervals = lagged_starts .- lagged_ends
    return durations[min_duration .< durations .< max_duration], intervals[min_interval .< intervals .< max_interval]
end

function extract_interval(timestamp_arr::Vector{Matrix{T}}; 
        flatten = true, kwargs...
    ) where T <: Real
    if flatten
        #In this case we don't necessarily need to preserve the structure data and can collapse all entries into one
        durations = T[]
        intervals = T[]
        for idx in 1:length(timestamp_arr)
            result = extract_interval(timestamp_arr[idx]; kwargs...)
            if !isnothing(result)
                push!(durations, result[1]...)
                push!(intervals, result[2]...)
            end
        end
        return durations, intervals
    else
        println("Not implemented")
    end
end
"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
    It takes in timestamps and returns the burst durations and the spikes per burst
A multiple dispatch of this function allows the max_interval to be calculated on a 3D array (x, y, and time) 
"""
function max_interval_algorithim(timestamps::Matrix{T}; 
        ISIstart::T = 500.0, ISIend::T = 500.0, IBImin::T = 1000.0, DURmin::T = 500.0, SPBmin::Int64 = 4, 
        verbose = false
    ) where T <: Real
    burst_timestamps = Tuple[]
    SPB_list = Float64[]
    if isempty(timestamps)
        if verbose >= 1
            println("No spikes detected")
        end
        return nothing
    else
        #Lets organize the spipkes into intervals spikes and not spikes
        results = extract_interval(timestamps)
        intervals = results[2]
        bursting = false
        burst_start_list = T[]
        burst_end_list = T[]
        burst_start = 0.0
        burst_end = 0.0
        SPB = 0
        idx = 1
        for i = 1:length(intervals)
            if bursting == false && intervals[i] <= ISIstart #If the cell does not start as bursting and the interval is under ISI start
                bursting = true #Begin the burst
                burst_start = timestamps[i, 1] #Record the burst start
            elseif bursting == true && intervals[i] >= ISIend || i == length(intervals) #If the cell is bursting, and the interval to the next spike is greater than ISI thresh
                bursting = false #The bursting can stop
                burst_end = timestamps[i, 2] #The burst end can be recorded
                if intervals[i] >= IBImin && (burst_end - burst_start) >= DURmin && SPB >= SPBmin #If the burst meets all the correct qualifications
                    if verbose
                        println("
                        Burst #$idx successfully added at timestamp : $burst_start -> $burst_end 
                            Duration: $(burst_end - burst_start) >  $DURmin  
                            Spikes per burst: $SPB > $SPBmin
                            IBI to burst #$(idx+1): $(intervals[i])
                            "
                            )
                    end
                    push!(burst_start_list, burst_start)
                    push!(burst_end_list, burst_end)
                    #push!(burst_timestamps, (burst_start, burst_end)) #Record it
                    push!(SPB_list, SPB)
                    SPB = 0
                    idx+=1
                elseif i == length(intervals) && (burst_end - burst_start) >= DURmin && SPB >= SPBmin
                    #a weird caveat, bursting has finished but interval has never cleared the ISIend
                    if verbose
                        println("
                        Burst #$idx successfully added at timestamp : $burst_start -> $burst_end 
                            Duration: $(burst_end - burst_start) >  $DURmin  
                            Spikes per burst: $SPB > $SPBmin
                            "
                            )
                    end
                    push!(burst_start_list, burst_start)
                    push!(burst_end_list, burst_end)
                    #push!(burst_timestamps, (burst_start, burst_end)) #Record it
                    push!(SPB_list, SPB)
                    SPB = 0
                    idx+=1
                else
                    if verbose
                        println("
                        Burst did not fit recommended qualities
                            Timestamp $idx: $burst_start -> $burst_end, 
                            DUR $idx: $(burst_end - burst_start) <  $DURmin 
                            SPB $idx: $SPB < $SPBmin
                            IBI $idx: $(intervals[i])
                            "
                        )
                    end                    
                end
            end
            if bursting == true
                SPB += 1
            end
        end
        

        if length(burst_start_list) > length(burst_end_list)
            #This algorithim usually leaves one last burst off because it has no end point. We can add this
            push!(burst_end_list, burst_start_list[end] + intervals[end])
        end        
        burst_timestamps = hcat(burst_start_list, burst_end_list)
        if isempty(burst_start_list)
            return nothing
        end
        return burst_timestamps, SPB_list
    end
end

function max_interval_algorithim(timestamp_arr::Vector{Matrix{T}}; kwargs...) where T <: Real
    #In this case we don't necessarily need to preserve the structure data and can collapse all entries into one
    bursts = Matrix{T}[]
    spd = T[]
    for idx in 1:length(timestamp_arr)
        
        result = max_interval_algorithim(timestamp_arr[idx]; kwargs...)
        if !isnothing(result)
            push!(bursts, result[1])
            push!(spd, result[2]...)
        end
    end
    return bursts, spd
end


function timeseries_analysis(sol::DiffEqBase.AbstractODESolution; 
        dt::Float64 = 100.0, Z::Int64 = 4,  
        max_spike_duration::Float64 = 10.0, max_burst_duration::Float64 = 10e5
    )
    thresholds = calculate_threshold(sol; Z = Z, dt = dt) #This takes really long
    spikes = get_timestamps(sol, thresholds)
    spike_durs, isi = extract_interval(spikes, max_duration = max_spike_duration)
    bursts, spb = max_interval_algorithim(spikes)
    burst_durs, ibi = extract_interval(bursts, max_duration = max_burst_duration)

    timestamps = Dict(
        "Spikes" => spikes,
        "Bursts" => bursts
    )

    data = Dict(
        "Thresholds" => thresholds,
        "SpikeDurs" => spike_durs, 
        "ISIs" => isi, 
        "BurstDurs" => burst_durs, 
        "IBIs" => ibi,
        "SpikesPerBurst" => spb
    )

    return timestamps, data
end


function timeseries_analysis(save_file::String, sol::DiffEqBase.AbstractODESolution;
        plot_histograms = true, kwargs...   
    )
    timestamps, data = timeseries_analysis(sol; kwargs...)
    bson("$(save_file)\\timestamps.bson", timestamps)
    bson("$(save_file)\\data.bson", data)

    if plot_histograms && !isempty(data["SpikeDurs"])
        sdur_hfit = fit(Histogram, data["SpikeDurs"], LinRange(0.0, 50.0, 1000))
        sdur_weights = sdur_hfit.weights/maximum(sdur_hfit.weights)
        sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]
        p1 = plot(sdur_edges, sdur_weights, xlabel = "Spike Duration (ms)")

        isi_hfit = fit(Histogram, data["ISIs"], LinRange(0.0, 100.0, 1000))
        isi_weights = isi_hfit.weights/maximum(isi_hfit.weights)
        isi_edges = collect(isi_hfit.edges[1])[1:length(isi_weights)]
        p2 = plot(isi_edges, isi_weights,  xlabel = "Spike Interval (s)", xformatter = x -> x/1000)

        bdur_hfit = fit(Histogram, data["BurstDurs"], LinRange(0.0, 2000.0, 1000))
        bdur_weights = bdur_hfit.weights/maximum(bdur_hfit.weights)
        bdur_edges = collect(bdur_hfit.edges[1])[1:length(bdur_weights)]
        p3 = plot(bdur_edges, bdur_weights, xlabel = "Burst Duration (s)",xformatter = x -> x/1000)

        ibi_hfit = fit(Histogram, data["IBIs"], LinRange(0.0, 60e3, 1000))
        ibi_weights = ibi_hfit.weights/maximum(ibi_hfit.weights)
        ibi_edges = collect(ibi_hfit.edges[1])[1:length(ibi_weights)]
        p4 = plot(ibi_edges, ibi_weights, xlabel = "Interburst Interval (s)", xformatter = x -> x/1000)
        
        p5 = histogram(data["Thresholds"], yaxis = :log, xlabel = "Voltage threshold")
        p6 = histogram(data["SpikesPerBurst"], yaxis = :log, xlabel = "Spiked per Burst")
        
        
        hist_plot = plot(
            p1, p2, p3, p4, p5, p6, 
            layout = grid(3,2), ylabel = "Counts", 
            legend = false
            )
        savefig(hist_plot, "$(save_file)\\histogram_plot.png")
    end

    BotNotify("{Wave} Finished running timeseries analysis")
    return timestamps, data
end

"""
Extract the waves
"""
function extract_waves(sol::DiffEqBase.AbstractODESolution, thresholds::Matrix{T};
        wave_min = 200 #This is the minimum number of points in the wave    
    ) where T<:Real
    nx = ny = round(Int64, sqrt(size(sol,1))) #This is the x and y dimension of the graph
    spike_array = (sol |> Array).>thresholds
    spike_array = reshape(spike_array, (nx, ny, size(sol,2)))
    markers = label_components(spike_array)
    #Maybe reshaping markers would help
    markers = reshape(markers, (nx*ny, size(markers,3)))
    n_markers = maximum(markers)
    good_markers = zeros(size(markers))
    for i = 1:n_markers #This simply prunes all unnecessary coordinates
        idxs = (markers .== i)
        if sum(idxs) >= wave_min #We can consider this a wave
            println(i)
            good_markers .+= idxs
        end
    end
    return Matrix{Bool}(good_markers.==1.0)
end



"""
For 3D arrays and functions, this will extract all of the bursts and convert it into a graphable array
"""
function extract_burstmap(spike_array::BitArray{3})
    burst_map = zeros(size(spike_array))
    burst_data = max_interval_algorithim(spike_array; dt = 1.0)
    for (x,y,data) in burst_data
        if !isnothing(data[1])
            for (rng_start, rng_stop) in data[1]
                burst_map[x, y, Int(rng_start):Int(rng_stop)] .= 1.0
            end
        end
    end
    return burst_map
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
