#### Time Scale Analysis ################################################################

"""
    calculate_threshold(sol::DiffEqBase.AbstractODESolution, Z::Int64)

Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation. 
If using a differential solution, make sure dt is set, otherwise the standard deviation will be unevenly sampled
"""
function calculate_threshold(vm_arr::AbstractArray; Z = 4, dims = -1)
    if dims == -1
        return [sum(vm_arr)/length(vm_arr) + Z*std(vm_arr)]
    else
        n = size(vm_arr, dims)
        mean = sum(vm_arr, dims = dims)./n
        dev = Z * std(vm_arr, dims = dims)
        return mean + dev #We want these all to come out as vectors vs matrices
    end
end

#This form of the function only calculates the analysis on a single variable
function calculate_threshold(sol::DiffEqBase.AbstractODESolution{T, N, S}, rng::Tuple{Float64, Float64};
        idxs::Union{Int64, AbstractArray{Int64}} = 1, Z=4.0, dt::T=100.0
    ) where {T, N, S}
    
    #println(N) #N Represents how many dimensions the data is 

    if N == 2 #This represents a single neuron simulation N = (variable, time)
        #We want to check how many dimensions the simulation is
        if idxs == -1 #This computes the analysis on all variables
            data_section = sol(rng[1]:dt:rng[2])
        else #This only looks at a single variable
            data_section = sol(rng[1]:dt:rng[2], idxs = idxs)
        end
        if length(size(data_section)) == 1 #The array is just a single solution (Var, Time)
            return calculate_threshold(data_section; Z = Z)
        else #This means the array has a n of (X, Time). We want to summarize by X
            return calculate_threshold(data_section; Z = Z, dims = 2)
        end
    elseif N == 4 #This represents a grid simulation N = (x, y, variable, time)
        trng = rng[1]:dt:rng[2]
        nx = size(sol, 1)
        ny = size(sol, 2)
        nt = length(trng)
        if isa(idxs, AbstractArray) #If you want to threshold multiple arrays
            data_section = zeros(nx, ny, length(idxs), nt)
            for (i, idx) in enumerate(idxs)
                #We want to pull out each index. Remember these are automatically flattened
                data_i = sol(trng, idxs=[1+(idx-1)*(nx*ny):(idx)*(nx*ny)...]) |> Array
                data_section[:, :, i, :] = reshape(data_i, (nx, ny, nt))
            end
        elseif idxs == -1
            data_section = sol(trng)
        else
            data_section = sol(trng, idxs = [1+(idxs-1)*(nx*ny):(idxs)*(nx*ny)...]) |> Array
            data_section = reshape(data_section, (nx, ny, 1, nt))
        end
        return calculate_threshold(data_section; Z=Z, dims=4)
    end
end

#This form just is a shortcut if you don't provide a time range
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
    diff_starts = findall(x -> x==-1, diff_vals) #This is a list of all the starting points in the array
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

#This dispatch calls the function if the time series is in the form of a range vs a vector
get_timestamps(spike_array, time_range::StepRangeLen{T, P, P}) where {T <: Real, P} = get_timestamps(spike_array, collect(time_range))

"""
If dt is set to Inf, the algorithim acts adaptive
"""
function get_timestamps(sol::DiffEqBase.AbstractODESolution{T, N, S}, threshold, rng::Tuple{T, T}; 
        idxs::Union{Int64, AbstractArray{Int64}} = 1,  
        dt::Union{T, Symbol} = :adaptive, 
        flatten::Bool = false,
        verbose::Bool = true
    ) where {T <: Real, N, S <: Vector} #When editing we need to focus on this section
    # Lets do this the integrated way
    # First we need to extract the spike array
    if N == 2 #This is for a single cell simulation
        nVar = size(sol, 1)
        nt = size(sol, 2)
        if dt == :adaptive #If the time range is adaptive use the dt of the solution
            time_rng = sol.t[rng[1].<sol.t.<rng[2]]
        else
            time_rng = rng[1]:dt:rng[2]
        end
    
        if idxs == -1
            @assert length(threshold) == nVar #We want to ensure the number of thresholds equals the number of variables
            timestamps = Vector{Matrix{T}}()
            for idx in 1:nVar
                data_select = sol(time_rng, idxs=idx) |> Array
                spike_array = Vector{Bool}(data_select .> threshold[idx]) #This always needs to be a Vector of Bools
                push!(timestamps, get_timestamps(spike_array, time_rng))
            end
            return timestamps
        elseif isa(idxs, AbstractArray)
            @assert length(threshold) == length(idxs) #We want to ensure the number of thresholds equals the number of variables
            timestamps = Vector{Matrix{T}}()
            for (i, idx) in enumerate(idxs)
                data_select = sol(time_rng, idxs=idx) |> Array
                spike_array = Vector{Bool}(data_select .> threshold[i]) #This always needs to be a Vector of Bools
                push!(timestamps, get_timestamps(spike_array, time_rng))
            end
            return timestamps
        else
            data_select = sol(time_rng, idxs=idxs) |> Array
            spike_array = Vector{Bool}(data_select .> threshold) #This always needs to be a Vector of Bools
            return get_timestamps(spike_array, time_rng)
        end
    elseif N == 4 #This is for a larger scale simulation
        if dt == :adaptive #If the time range is adaptive use the dt of the solution
            time_rng = sol.t[rng[1].<sol.t.<rng[2]]
        else
            time_rng = rng[1]:dt:rng[2]
        end
        nx = size(sol, 1)
        ny = size(sol, 2)
        nVar = size(sol, 3)
        nt = length(time_rng)
        @assert nx == size(threshold, 1)
        @assert ny == size(threshold, 2)

        #We should collapse all thresholds and flatten all variables
        if idxs == -1 #This extracts all the variables
            @assert nVar == size(threshold, 3)
            thresh = reshape(thresholds, nx, ny, nVar)
            data_select = sol(time_rng) |> Array
            println(size(data_select))
            timestamps = Array{Matrix{T}}()
            #for idx in 1:nVar
            #    println(idx)
            #end
        elseif isa(idxs, AbstractArray)
            data_select = zeros(nx, ny, length(idxs), nt)
            for (i, idx) in enumerate(idxs)
                data_i = sol(time_rng, idxs=[1+(idx-1)*(nx*ny):(idx)*(nx*ny)...]) |> Array
                data_select[:, :, i, :] = reshape(data_i, nx, ny, nt)
            end
            println(size(data_select))
        else
            thresh = reshape(threshold, nx*ny)
            data_select = sol(time_rng, idxs=[1+(idxs-1)*(nx*ny):(idxs)*(nx*ny)...]) |> Array
            timestamps = Vector{Matrix{T}}(undef, nx*ny)
            for i in 1:size(data_select, 1)
                if verbose 
                    println("Analyzing cell $i of $(size(data_select,1))")
                end
                spike_array = Vector{Bool}(data_select[i, :] .> thresh[i]) #This always needs to be a Vector of Bools
                stamps = get_timestamps(spike_array, time_rng)
                timestamps[i] = stamps
            end
        end
        return timestamps
    end
end

get_timestamps(sol::DiffEqBase.AbstractODESolution{T, N, S}, threshold; kwargs...) where {T <: Real, N, S <: Vector} = get_timestamps(sol, threshold, (sol.t[1], sol.t[end]); kwargs...)

function get_timestamps(sol::DiffEqBase.AbstractODESolution{T,N,S}; idxs = 1, Z::Int64 = 4, kwargs...) where {T <: Real, N, S <: Vector} 
    threshold = calculate_threshold(sol; Z = Z, idxs = idxs)
    return get_timestamps(sol, threshold, idxs = idxs, kwargs...)
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

function extract_interval(timestamp_arr::Vector{Union{Matrix{T}, Nothing}}; 
        flatten = true, kwargs...
    ) where T <: Real
    #In this case we don't necessarily need to preserve the structure data and can collapse all entries into one
    durations = T[]
    intervals = T[]
    for (idx, tstamps) in enumerate(timestamp_arr)
        if !isnothing(tstamps)
            result = extract_interval(tstamps; kwargs...)
            if !isnothing(result)
                push!(durations, result[1]...)
                push!(intervals, result[2]...)
            end
        end
    end
    return durations, intervals
end

#This might be an error, but this function doesn't seem to work unless you have two seperate versions
function extract_interval(timestamp_arr::Vector{Matrix{T}}; 
        flatten = true, kwargs...
    ) where T <: Real
    #In this case we don't necessarily need to preserve the structure data and can collapse all entries into one
    durations = T[]
    intervals = T[]
    for (idx, tstamps) in enumerate(timestamp_arr)
        result = extract_interval(tstamps; kwargs...)
        if !isnothing(result)
            push!(durations, result[1]...)
            push!(intervals, result[2]...)
        end
    end
    return durations, intervals
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
    It takes in timestamps and returns the burst durations and the spikes per burst
"""
function max_interval_algorithim(timestamps::Matrix{T}; 
        ISIstart::T = 500.0, ISIend::T = 500.0, IBImin::T = 1000.0, DURmin::T = 50.0, SPBmin::Int64 = 4, 
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
    bursts = Union{Matrix{T}, Nothing}[]
    spd = Union{T, Nothing}[]
    for idx in 1:length(timestamp_arr)
        if isassigned(timestamp_arr, idx)
            result = max_interval_algorithim(timestamp_arr[idx]; kwargs...)
            if !isnothing(result)
                push!(bursts, result[1])
                push!(spd, result[2]...)
            else
                push!(bursts, nothing)
                push!(spd, nothing)
            end
        end
    end
    return bursts, spd
end

function timeseries_analysis(t, vm_array;
    timestamps_only = false, Z::Int64=4, 
    max_spike_duration::Float64=50.0, max_spike_interval = 100,
    max_burst_duration::Float64=10e5, max_burst_interval = 10e5,
    verbose=false
) #where {T, N}
    N = length(size(vm_array))
    #println(N)
    if verbose
        print("[$(now())]: Extracting the thresholds... ")
    end
    if N==1
        thresholds = calculate_threshold(vm_array)
        spike_array = Vector{Bool}(vm_array .> thresholds)
        #return thresholds
    elseif N == 2
        thresholds = calculate_threshold(vm_array, dims=2)
        spike_array = Matrix{Bool}(vm_array .> thresholds)
    end
    #println(spike_array |> typeof)
    if verbose
        println("Completed")
    end
    spikes = get_timestamps(spike_array, t)
    res = max_interval_algorithim(spikes)
    #println(res)
    #println(isnothing(res))
    if isnothing(res)
        bursts = spb = nothing
    else
        bursts, spb = res
    end

    timestamps = Dict(
        "Spikes" => spikes,
        "Bursts" => bursts
    )
    if timestamps_only
        return timestamps
    else
        if verbose
            print("[$(now())]: Extracting the Data... ")
        end
        burst_durs, ibi = extract_interval(bursts, max_duration=max_burst_duration, max_interval=max_burst_interval)
        spike_durs, isi = extract_interval(spikes, max_duration=max_spike_duration, max_interval=max_spike_interval)
        data = Dict(
            "Time" => t,
            "DataArray" => vm_array,
            "Thresholds" => thresholds,
            "SpikeDurs" => spike_durs,
            "ISIs" => isi,
            "BurstDurs" => burst_durs,
            "IBIs" => ibi,
            "SpikesPerBurst" => spb
        )
        if verbose
            println("Complete")
        end
        return timestamps, data
    end
end

function timeseries_analysis(t::AbstractArray{T}, vm_array::Array{T, N}, save_file::String;
    tstamps_name="timestamps", data_name="data",
    verbose = false, 
    kwargs...
) where {T, N}
    timestamps, data = timeseries_analysis(t, vm_array; kwargs...)
    if verbose
        print("[$(now())]: Saving data... ")
    end
    #Uncomment to use BSON file format
    #bson("$(save_file)\\timestamps.bson", timestamps)
    #bson("$(save_file)\\data.bson", data)
    #Uncomment to use JLD2 to save the packages
    save("$(save_file)/$(tstamps_name).jld2", timestamps)
    save("$(save_file)/$(data_name).jld2", data)
    if verbose
        println("Complete")
    end
    return timestamps, data
end


function timeseries_analysis(sol::DiffEqBase.AbstractODESolution{T,N,S};
    dt=1.0, kwargs...
) where{T, N, S}
    t = sol.t[1]:dt:sol.t[end]
    vm_array = sol(t) |> Array
    return timeseries_analysis(t, vm_array; kwargs...)
end

function timeseries_analysis(sol::DiffEqBase.AbstractODESolution, save_file::String; 
        tstamps_name =  "timestamps", data_name = "data", 
        kwargs...
    )
    timestamps, data = timeseries_analysis(sol; kwargs...)
    print("[$(now())]: Saving data... ")
    #Uncomment to use BSON file format
    #bson("$(save_file)\\timestamps.bson", timestamps)
    #bson("$(save_file)\\data.bson", data)
    #Uncomment to use JLD2 to save the packages
    save("$(save_file)/$(tstamps_name).jld2", timestamps)
    save("$(save_file)/$(data_name).jld2", data)
    println("Complete")
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