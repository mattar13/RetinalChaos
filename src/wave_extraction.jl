#To extend some utilites you need to explicitly import them
import ePhys: calculate_threshold, get_timestamps, timeseries_analysis, max_interval_algorithim

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

function timeseries_analysis(sol::DiffEqBase.AbstractODESolution{T,N,S};
    dt=1.0, idxs = -1, kwargs...
) where{T, N, S}
    t = sol.t[1]:dt:sol.t[end]
    if idxs == -1
        vm_array = sol(t) |> Array
    else
        vm_array = sol(t, idxs = idxs) |> Array
    end
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

#function extract_data()


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
function extract_burstmap(nx, ny, nt, burst_timestamps)
    burst_map = falses(nx*ny,nt)
    for (cell_idx, cell) in enumerate(burst_timestamps)
        #println(cell_idx)
        #println(cell)
        if !isnothing(cell)
            #println(cell)
            for idx in 1:size(cell,1)
                rng_start, rng_stop = cell[idx, :]
                if rng_start == 0
                    rng_start = 1
                end
                burst_map[cell_idx, Int(rng_start):Int(rng_stop)] .= true
            end
        end
    end
    burst_map = reshape(burst_map, nx, ny, nt)
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

function WaveSegmentation(spikes; min_wavesize = 20.0)
    #nx = ny = round(Int64, sqrt(size(spikes, 1)))
    tspan = size(spikes, 2)
    spike_map = spikes#reshape(spikes, (nx, ny, tspan)) #Reshape the spike array to be 3D
    wave_segs = label_components(spike_map) #Label all of the 3D components
    n_waves = maximum(wave_segs) #Pull out the maximum number of waves
    is_wave = falses(n_waves)
    wave_areas = zeros(n_waves) #Record all of the wave sizes 
    wave_times = zeros(n_waves)
    wave_distances = zeros(n_waves)
    wave_velocities = zeros(n_waves)
    for n = 1:n_waves #For each labeled wave
        print("Measuring wave $n of $n_waves: ")
        pᵢ = findall(wave_segs.==n) #Return a bitarray of all the points
        wave_areas[n] = length(pᵢ)
        println("Wave size is $(wave_areas[n])")
        #return pᵢ
        (xs, ys, zs) = RetinalChaos.extract_points(pᵢ)
        if (maximum(xs) - minimum(xs)) * (maximum(ys) - minimum(ys)) == 0.0
            #This means that the resulting wave is a burst and doesn't spread
            println(" Resulting wave is only a burst")
        else
            println(" Resulting wave is a wave")
            is_wave[n] = true
            wave_times[n] = maximum(zs) - minimum(zs)
            wave_distances[n] = (maximum(xs) - minimum(xs)) * (maximum(ys) - minimum(ys))
            wave_velocities[n] = wave_times[n] * wave_distances[n]
        end
    end
    df = DataFrame(
        WaveN = 1:n_waves, IsWave = is_wave,
        WaveArea = wave_areas, WaveTime = wave_times, 
        WaveDistances = wave_distances, WaveVelocities = wave_velocities
    )
    return df
end

function WaveSegmentation(spikes, save_loc::String; kwargs...)
    df = WaveSegmentation(spikes; kwargs...)
    only_bursts = df |> @filter(_.IsWave == true) |> DataFrame

    save_path = "$save_loc\\wave_data.xlsx"
    if isfile(save_path)
        println("File already exists")
        rm(save_path)
        println("Old file deleted")
    end
    XLSX.writetable(save_path, "WaveStats" => df, "Waves" => only_bursts)
    println("New file written")
    return df
end