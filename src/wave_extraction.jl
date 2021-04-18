#### Time Scale Analysis ################################################################
using DiffEqBase
"""
Calculate threshold
Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation 
"""
calculate_threshold(vm_arr::AbstractArray; Z::Int64 = 4) = sum(vm_arr)/length(vm_arr) + Z*std(vm_arr)

function calculate_threshold(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T, T}; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T= 0.1,
    ) where T <: Real
    data_section = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    return calculate_threshold(data_section; Z = Z)
end

function calculate_threshold(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real 
    data_section = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    calculate_threshold(data_section; Z = Z)
end

"""
This function acts to calculate the distance between points in a single BitArray. 
Very rarely is the first point part of a spike (in which case there is a fallback), 
and because of this clip is set to remove the first interval. 
"""
function count_intervals(spike_trace::BitArray{1}; 
        clip = 2, clip_end = 0
    )
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
    isi[clip:end-clip_end]
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
function get_timestamps(spike_array::BitArray{1}; 
        dt = 0.1, verbose = 0
    )
    idx_array = findall(x -> x==1, spike_array)
    points = Tuple[]
    if !isempty(idx_array)
        start_point = idx_array[1]
        end_point = idx_array[2]
        for i in 1:length(idx_array)-1
            if (idx_array[i+1] - idx_array[i]) != 1
                end_point = idx_array[i]
                push!(points, (start_point*dt, end_point*dt))
                start_point = idx_array[i+1]
            end
        end
    end
    points
end

function get_timestamps(spike_array::BitArray{2}; 
        dt = 0.1
    )
    nx, tsteps = size(spike_array)
    timestamps = Tuple[]
    for x = 1:nx
        stamps = get_timestamps(spike_array[x,:]; dt = dt)
        push!(timestamps, stamps)
    end
    timestamps
end

function get_timestamps(spike_array::BitArray{3}; 
        dt = 0.1
    )
    nx, ny, tsteps = size(spike_array)
    timestamps = Tuple[]
    for x = 1:nx
        for y = 1:ny
            stamps = get_timestamps(spike_array[x,y,:]; dt = dt)
            push!(timestamps, stamps)
        end
    end
    timestamps
end

function get_timestamps(sol::DiffEqBase.AbstractODESolution, threshold::T, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1
    ) where T <: Real
    #First we need to extract the spike array
    data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    get_timestamps(spike_array; dt = dt)
end

# For if the threshold has not been calculated
function get_timestamps(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T,T}; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real
    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    get_timestamps(spike_array; dt = dt)
end

# For if no range has been provided but a threshold has
function get_timestamps(sol::DiffEqBase.AbstractODESolution, threshold::T; 
        idx::Int64 = 1, dt::T = 0.1
    ) where T <: Real
    #First we need to extract the spike array
    data_select = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    get_timestamps(spike_array; dt = dt)
end

# For if no range has been provided
function get_timestamps(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real
    threshold = calculate_threshold(sol; idx = idx, Z = Z, dt = dt)
    data_select = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    get_timestamps(spike_array; dt = dt)
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
A multiple dispatch of this function allows the max_interval to be calculated on a 3D array (x, y, and time) 
"""
function max_interval_algorithim(spike_array::BitArray{1}; 
        ISIstart::Int64 = 500, ISIend::Int64 = 500, IBImin::Int64 = 1000, DURmin::Int64 = 100, SPBmin::Int64 = 4, 
        dt::T = 0.1, verbose = false
    ) where T <: Real
    burst_timestamps = Array{Tuple,1}([])
    DUR_list = Array{Float64,1}([])
    SPB_list = Array{Float64,1}([])
    IBI_list = Array{Float64,1}([])
    #Detect the spikes first
    timestamps = get_timestamps(spike_array; dt = dt) #Add in arguments for sweeps later
    if isempty(timestamps)
        if verbose >= 1
            println("No spikes detected")
        end
        return fill(nothing, 4)
    else
        #println("Times detected")
        #Lets organize the spipkes into intervals spikes and not spikes
        spike_durs = map(i -> timestamps[i][2]-timestamps[i][1], 1:length(timestamps))
        intervals = map(i -> timestamps[i][1] - timestamps[i-1][2], 2:length(timestamps))
        #intervals = count_intervals(spike_array) .* dt
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
        if DUR >= DURmin && SPB >= SPBmin && bursting == true
            if verbose
                println("Timestamp  $idx: $burst_start -> $(timestamps[end][2]), DUR $idx: $DUR, SPB $idx: $SPB, IBI $idx: Unknown")
            end
            push!(burst_timestamps, (burst_start, timestamps[end][2]))
            push!(DUR_list, DUR)
            push!(SPB_list, SPB)
        end
        return burst_timestamps, DUR_list, SPB_list, IBI_list
    end
end

function max_interval_algorithim(spike_array::BitArray{3}; kwargs...)
    nx, ny, tsteps = size(spike_array)
    data_array = Tuple[]
    for x = 1:nx
        for y = 1:ny
            data = max_interval_algorithim(spike_array[x,y,:]; kwargs...)
            push!(data_array, (x, y, data))
        end
    end
    data_array
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution, threshold::T, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1,
        kwargs...
    ) where T <: Real

    data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    return max_interval_algorithim(spike_array; dt = dt, kwargs...)
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T,T};         
        idx::Int64 = 1, dt::Float64 = 0.1, Z::Int64 = Z,
        kwargs...
    ) where T <: Real

    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    return max_interval_algorithim(spike_array; dt = dt, kwargs...)
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution, threshold::T; 
        idx::Int64 = 1, dt::T = 0.1, 
        kwargs...
    ) where T <: Real
    #First we need to extract the spike array
    data_select = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    return max_interval_algorithim(spike_array; dt = dt, kwargs...)
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, dt::T = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real
    #First we need to extract the spike array
    threshold = calculate_threshold(sol; idx = idx, Z = Z, dt = dt)
    data_select = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    return max_interval_algorithim(spike_array; dt = dt, kwargs...)
end

"""
Timescale analysis

returns Spike Durations, Burst Durations, Interburst Intervals
"""
#Basis for all other timescale analysis
function timescale_analysis(spike_array::BitArray{1}; 
        dt::T = 0.1,
        kwargs...
    ) where T <: Real

    if !any(spike_array .== 1.0)
        if verbose >= 1
            println("No spikes detected")
        end
        return fill(NaN, 3)
    else
        timestamps = get_timestamps(spike_array; dt = dt);
        durations = map(x -> x[2]-x[1], timestamps)
        if durations == Any[]
            #This essentially means that no spikes are detected, therefore no bursts occur
            return fill(NaN, 3)
        end
        burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array; dt = dt, kwargs...);
        return durations, dur_list, ibi_list
    end
end

function timescale_analysis(sol::DiffEqBase.AbstractODESolution, threshold::T, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1,
        kwargs...
    ) where T <: Real
    data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    timescale_analysis(spike_array; dt = dt, kwargs...)
end

function timescale_analysis(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real

    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    timescale_analysis(spike_array; dt = dt, kwargs...)
end

function timescale_analysis(sol::DiffEqBase.AbstractODESolution, threshold::T; 
        idx::Int64 = 1, dt::T = 0.1, 
        kwargs...
    ) where T <: Real
    #First we need to extract the spike array
    data_select = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    timescale_analysis(spike_array; dt = dt, kwargs...)
end

function timescale_analysis(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, dt::T = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real
    #First we need to extract the spike array
    threshold = calculate_threshold(sol; idx = idx, Z = Z, dt = dt)
    data_select = sol(sol.t[1]:dt:sol.t[end], idxs = idx) |> Array
    spike_array = (data_select .> threshold)
    timescale_analysis(spike_array; dt = dt, kwargs...)
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
