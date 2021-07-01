#### Time Scale Analysis ################################################################

"""
    calculate_threshold(sol::DiffEqBase.AbstractODESolution, Z::Int64)

Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation 
"""
calculate_threshold(vm_arr::AbstractArray; Z::Int64 = 4) = [sum(vm_arr)/length(vm_arr) + Z*std(vm_arr)]

function calculate_threshold(sol::DiffEqBase.AbstractODESolution{T, N, S}, rng::Tuple{Float64, Float64}; 
        idx::Int64 = 1, Z::Int64 = 4, dt::Float64 = 0.1,
    ) where {T, N, S}
    # We need to convert the dt into the correct form
    println("here")
    dt = convert(T, dt)
    #We want to check how many dimensions the simulation is 
    if length(size(sol.prob.u0)) == 1
        data_section = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
        #println(data_section |> size)
        return calculate_threshold(data_section; Z = Z)
    else  
        n = length(collect(rng[1]:dt:rng[2]))
        mean = sum(sol(rng[1]:dt:rng[2]), dims = 2)./n
        println("GPU solution")
        dev = Z * std(sol(rng[1]:dt:rng[2]))
        return mean .+ dev
    end
end

function calculate_threshold(sol::DiffEqBase.AbstractODESolution{T, N, S}; kwargs...) where {T, N, S} 
    t_rng = (sol.t[1]|> Float64, sol.t[end] |> Float64)
    println(t_rng)
    ans = calculate_threshold(sol, t_rng; kwargs...)
    
    return ans
end

"""
This function returns all the time stamps in a spike or burst array
"""
function get_timestamps(spike_array::BitArray{1}; 
        dt = 0.1, verbose = 0
    )
    idx_array = findall(x -> x==1, spike_array)
    points = Tuple[]
    if length(idx_array) > 1
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

function get_timestamps(sol::DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1
    ) where T <: Real
    # Lets do this the integrated way
    #First we need to extract the spike array
    if length(size(sol.prob.u0)) == 1
        data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
        spike_array = (data_select .> threshold[1])
        get_timestamps(spike_array; dt = dt)
    else
        #For each data trace we have a corresponding threshold
        timestamps = Tuple[]
        @showprogress "Extracting timestamps: " for c_idx in 1:size(sol,1)
            data_select = sol(rng[1]:dt:rng[2], idxs = c_idx) |> Array
            spike_array = (data_select .> threshold[c_idx])
            stamps = get_timestamps(spike_array; dt = dt)
            push!(timestamps, stamps...)
        end
        return timestamps
    end
    
end

# For if the threshold has not been calculated
function get_timestamps(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T,T}; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real
    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    get_timestamps(sol, threshold, rng; dt = dt)
end

# For if no range has been provided but a threshold has
function get_timestamps(sol::DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}; 
        idx::Int64 = 1, dt::T = 0.1
    ) where T <: Real
    get_timestamps(sol, threshold, (sol.t[1], sol.t[end]); idx = idx, dt = dt)
end

# For if no range has been provided
function get_timestamps(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, Z::Int64 = 4, dt::T = 0.1
    ) where T <: Real
    threshold = calculate_threshold(sol; idx = idx, Z = Z, dt = dt)
    get_timestamps(sol, threshold, (sol.t[1], sol.t[end]); idx = idx, dt = dt)
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
A multiple dispatch of this function allows the max_interval to be calculated on a 3D array (x, y, and time) 
"""
function max_interval_algorithim(spike_array::BitArray{1}; 
        ISIstart::Int64 = 500, ISIend::Int64 = 500, IBImin::Int64 = 1000, DURmin::Int64 = 100, SPBmin::Int64 = 4, 
        dt::T = 0.1, verbose = false
    ) where T <: Real
    burst_timestamps = Tuple[]
    DUR_list = Float64[]
    SPB_list = Float64[]
    IBI_list = Float64[]
    #Detect the spikes first
    timestamps = get_timestamps(spike_array; dt = dt) #Add in arguments for sweeps later
    if isempty(timestamps)
        if verbose >= 1
            println("No spikes detected")
        end
        return nothing
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

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1,
        kwargs...
    ) where T <: Real
    if length(size(sol.prob.u0)) == 1
        data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
        spike_array = (data_select .> threshold[1])
        return max_interval_algorithim(spike_array; dt = dt, kwargs...)
    else
        burst_timestamps = Tuple[]
        DUR_list = Float64[]
        SPB_list = Float64[]
        IBI_list = Float64[]
        @showprogress "Running max_interval_algorithim: " for c_idx in 1:size(sol,1)
            data_select = sol(rng[1]:dt:rng[2], idxs = c_idx) |> Array
            spike_array = (data_select .> threshold[c_idx])
            MI_data = max_interval_algorithim(spike_array; kwargs...)
            if !isnothing(MI_data)
                push!(burst_timestamps, MI_data[1]...)
                push!(DUR_list, MI_data[2]...)
                push!(SPB_list, MI_data[3]...)
                push!(IBI_list, MI_data[4]...)
            end
        end
        return burst_timestamps, DUR_list, SPB_list, IBI_list
    end
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution, rng::Tuple{T,T};         
        idx::Int64 = 1, dt::Float64 = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real

    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    return max_interval_algorithim(sol, threshold, rng; dt = dt, kwargs...)
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}; 
        idx::Int64 = 1, dt::T = 0.1, 
        kwargs...
    ) where T <: Real
    return max_interval_algorithim(sol, threshold, (sol.t[1]. sol.t[end]); dt = dt, kwargs...)
end

function max_interval_algorithim(sol::DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, dt::T = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real
    #First we need to extract the spike array
    threshold = calculate_threshold(sol; idx = idx, Z = Z, dt = dt)
    return max_interval_algorithim(sol, threshold, (sol.t[1], sol.t[end]); dt = dt, kwargs...)
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
        #println("No spikes detected")
        return nothing
    else
        timestamps = get_timestamps(spike_array; dt = dt);
        spike_durations = map(x -> x[2]-x[1], timestamps)
        if spike_durations == Any[]
            #This essentially means that no spikes are detected, therefore no bursts occur
            return nothing
        end
        burst_idxs, burst_durations, spb_list, ibis = max_interval_algorithim(spike_array; dt = dt, kwargs...);
        return spike_durations, burst_durations, ibis
    end
end

function timescale_analysis(sol::DiffEqBase.DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1,
        kwargs...
    ) where T <: Real
    if length(size(sol.prob.u0)) == 1
        data_select = sol(rng[1]:dt:rng[2], idxs = idx) |> Array
        spike_array = (data_select .> threshold[1])
        return timescale_analysis(spike_array; dt = dt, kwargs...)
    else
        spike_durations = Float64[]
        burst_durations = Float64[]
        ibis = Float64[]
        @showprogress "Conducting timescale analysis: " for c_idx in 1:size(sol,1)
            data_select = sol(rng[1]:dt:rng[2], idxs = c_idx) |> Array
            spike_array = (data_select .> threshold[c_idx])
            ts_data = timescale_analysis(spike_array; dt = dt, kwargs...)
            if !isnothing(ts_data)
                push!(spike_durations, ts_data[1]...)
                push!(burst_durations, ts_data[2]...)
                push!(ibis, ts_data[3]...)
            end
        end
        return spike_durations, burst_durations, ibis
    end
    
end

function timescale_analysis(sol::DiffEqBase.DiffEqBase.AbstractODESolution, rng::Tuple{T,T}; 
        idx::Int64 = 1, dt::Float64 = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real

    threshold = calculate_threshold(sol, rng; idx = idx, Z = Z, dt = dt)
    timescale_analysis(sol, threshold, rng; dt = dt, kwargs...)
end

function timescale_analysis(sol::DiffEqBase.DiffEqBase.AbstractODESolution, threshold::AbstractArray{T}; 
        idx::Int64 = 1, dt::T = 0.1, 
        kwargs...
    ) where T <: Real
    timescale_analysis(sol, threshold, (sol.t[1], sol.t[end]); dt = dt, kwargs...)
end

function timescale_analysis(sol::DiffEqBase.DiffEqBase.AbstractODESolution; 
        idx::Int64 = 1, dt::T = 0.1, Z::Int64 = 4,
        kwargs...
    ) where T <: Real
    #First we need to extract the spike array
    threshold = calculate_threshold(sol; idx = idx, Z = Z, dt = dt)
    timescale_analysis(sol, threshold, (sol.t[1], sol.t[end]); dt = dt, kwargs...)
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
