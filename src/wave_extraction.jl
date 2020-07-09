#### Time Scale Analysis ################################################################

"""
Calculate threshold
Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation 
"""
calculate_threshold(vm_arr::AbstractArray where T; Z::Int64 = 4) = sum(vm_arr)/length(vm_arr) + Z*std(vm_arr)

"""
This dispatch calculates the threshold on a JLD file
"""
function calculate_threshold(filename::String; Z::Int64 = 4)
    tstamps = jldopen("test.jld", "r") do file
        read(file, "time")
    end
    n_points = nx*ny*length(tstamps)
    avg = 0.0
    covar = 0.0
    thresh = jldopen(filename, "r") do file
        for t in tstamps
            if t%10000==0.0
                #println(t)
            end
            if t == 0.0
                arr = read(file, "$(t)")[:,:,1]
                avg += sum(arr)
            else
                arr = read(file, "$(t)")
                avg += sum(arr)
            end
        end
        #Calculate the average
        avg /= n_points
        println(avg)
        for t in tstamps
            if t%10000==0.0
                #println(t)
            end
            if t == 0.0
                arr = read(file, "$(t)")[:,:,1]
                covar += sum((arr .- avg).^2)
            else
                arr = read(file, "$(t)")
                covar += sum((arr .- avg).^2)
            end
        end
        #Calculate the standard deviation
        std = sqrt(covar/n_points-1)
        #Calculate the threshold
        avg + Z*std
    end
    thresh
end

"""
This function acts to calculate the distance between points in a single BitArray. 
Very rarely is the first point part of a spike (in which case there is a fallback), 
and because of this clip is set to remove the first interval. 
"""
function count_intervals(spike_trace::BitArray{1}; clip = 2, clip_end = 0)
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
function get_timestamps(spike_array::BitArray{1}; dt = 1.0, verbose = 0)
    intervals = count_intervals(spike_array) .* dt
    durations = count_intervals(spike_array .!= 1.0) .* dt
    first_point = (findfirst(x -> x==1, spike_array)-1|>Float64) * dt
    current_point = first_point
    #If we land on the rare occasion where the last point of the spike array is true, we will have an incomplete final interval and interval and duration will be the same. 
    if length(intervals) == length(durations)
        points = Tuple[]
        for idx = 1:length(intervals)
            current_point += durations[idx] 
            push!(points, (current_point, current_point + durations[idx]))
            current_point += intervals[idx] 
        end
        return points        
    else
        points = Tuple[(current_point, current_point+durations[1])]
        current_point += durations[1]
        for idx = 1:length(intervals)
            current_point += intervals[idx] 
            push!(points, (current_point, current_point + durations[idx+1]))
            current_point += durations[idx+1] 
        end
        return points
    end 
end

function get_timestamps(spike_array::BitArray{3}; dt = 1.0)
    nx, ny, tsteps = size(spike_array)
    timestamps = Tuple[]
    for x = 1:nx
        for y = 1:ny
            stamps = get_timestamps(spike_array[x,y,:]; dt = dt)
            push!(timestamps, (x, y, stamps))
        end
    end
    timestamps
end

"""
This function uses the Maximum Interval Sorting method to sort bursts in a single trace. 
A multiple dispatch of this function allows the max_interval to be calculated on a 3D array (x, y, and time) 
"""
function max_interval_algorithim(spike_array::BitArray{1}; ISIstart = 500, ISIend = 500, IBImin = 1000, DURmin = 100, SPBmin = 4, dt = 1.0, verbose = false)
    burst_timestamps = Array{Tuple,1}([])
    DUR_list = Array{Float64,1}([])
    SPB_list = Array{Float64,1}([])
    IBI_list = Array{Float64,1}([])
    if !any(spike_array.==1.0)
        if verbose >= 1
            println("No spikes detected")
        end
        return fill(nothing, 4)
    else
        intervals = count_intervals(spike_array) .* dt
        timestamps = get_timestamps(spike_array; dt = dt)
        if length(timestamps) == 0
            #Somehow there is a weird case where there are spikes, but none that matter"
            return fill(nothing, 4)
        end
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

function max_interval_algorithim(spike_array::BitArray{3}; dt = 1.0)
    nx, ny, tsteps = size(spike_array)
    data_array = Tuple[]
    for x = 1:nx
        for y = 1:ny
            data = max_interval_algorithim(spike_array[x,y,:]; dt = dt)
            push!(data_array, (x, y, data))
        end
    end
    data_array
end

function timescale_analysis(vm_trace::Array{Float64,1}; dt = 10.0, verbose = 0, mode = 1)
    sim_thresh = calculate_threshold(vm_trace)
    spike_array = (vm_trace .> sim_thresh);
    if !any(spike_array .== 1.0)
        if verbose >= 1
            println("No spikes detected")
        end
        return fill(NaN, 6)
    else
        timestamps = get_timestamps(spike_array; dt = dt);
        durations = map(x -> x[2]-x[1], timestamps)
        burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array; verbose = (verbose>=2), dt = dt);
        
        avg_spike_dur = sum(durations)/length(durations)
        std_spike_dur = std(durations)
        avg_burst_dur = sum(dur_list)/length(dur_list);
        std_burst_dur = std(dur_list);
        avg_ibi = sum(ibi_list)/length(ibi_list)
        std_ibi = std(ibi_list)      

        if verbose >= 1
            println("The average spike duration is $(round(avg_spike_dur, digits = 2)) ms +- $(round(std_spike_dur, digits = 2)) ms")
            println("The average burst duration is $(round(avg_burst_dur;digits = 2)) ms +- $(round(std_burst_dur;digits = 2)) ms")
            println("The average interburst interval is $(round(avg_ibi;digits = 2)) ms +- $(round(std_ibi;digits = 2)) ms")
        end
        
        if mode == 1
            return [avg_spike_dur, std_spike_dur, avg_burst_dur, std_burst_dur, avg_ibi, std_ibi]
        elseif mode == 2
            #Mode 2 returns lists of spike durations, bursts and ibis
            return durations, dur_list, ibi_list
        end

    end
end

function timescale_analysis(vm_trace::Array{Float64,3}; dt = 10.0, verbose = 0, mode = 1)
    nx, ny, tsteps = size(vm_trace);
    thresh = calculate_threshold(vm_trace)
    spike_array = vm_trace .> thresh;
    spike_durations = Float64[]
    burst_durations = Float64[]
    IBIs = Float64[]
    for x = 1:nx
        for y = 1:ny
            if !any(spike_array[x,y,:] .== 1.0)
                if verbose >= 1
                    println("No spikes detected")
                end
            else
                timestamps = get_timestamps(spike_array[x,y,:]; dt = dt);
                durations = map(x -> x[2]-x[1], timestamps);
                push!(spike_durations, durations...)
                burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array[x,y,:]; verbose = (verbose>=2), dt = dt);
                if dur_list != nothing
                    push!(burst_durations, dur_list...)
                    push!(IBIs, ibi_list...)
                end
            end
        end
    end
    avg_spike_dur = sum(spike_durations)/length(spike_durations)
    std_spike_dur = std(spike_durations)
    avg_burst_dur = sum(burst_durations)/length(burst_durations);
    std_burst_dur = std(burst_durations);
    avg_ibi = sum(IBIs)/length(IBIs)
    std_ibi = std(IBIs)
    if verbose >= 1
        println("The average spike duration is $(round(avg_spike_dur, digits = 2)) ms +- $(round(std_spike_dur, digits = 2)) ms")
        println("The average burst duration is $(round(avg_burst_dur;digits = 2)) ms +- $(round(std_burst_dur;digits = 2)) ms")
        println("The average interburst interval is $(round(avg_ibi;digits = 2)) ms +- $(round(std_ibi;digits = 2)) ms")
    end
    if mode == 1
        return [avg_spike_dur, std_spike_dur, avg_burst_dur, std_burst_dur, avg_ibi, std_ibi]
    elseif mode == 2
        #Mode 2 returns lists of spike durations, bursts and ibis
        return spike_durations, burst_durations, IBIs
    end
end

#This dispatch calculates the timescale analysis on a already thresholded array
function timescale_analysis(spike_array::BitArray{3}; dt = 1.0, verbose = 0, mode = 1)
    nx, ny, tsteps = size(spike_array);
    spike_durations = Float64[]
    burst_durations = Float64[]
    IBIs = Float64[]
    for x = 1:nx
        for y = 1:ny
            if !any(spike_array[x,y,:] .== 1.0)
                if verbose >= 1
                    println("No spikes detected")
                end
            else
                timestamps = get_timestamps(spike_array[x,y,:]; dt = dt);
                durations = map(x -> x[2]-x[1], timestamps);
                push!(spike_durations, durations...)
                burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array[x,y,:]; verbose = (verbose>=2), dt = dt);
                if dur_list != nothing
                    push!(burst_durations, dur_list...)
                    push!(IBIs, ibi_list...)
                end
            end
        end
    end
    avg_spike_dur = sum(spike_durations)/length(spike_durations)
    std_spike_dur = std(spike_durations)
    avg_burst_dur = sum(burst_durations)/length(burst_durations);
    std_burst_dur = std(burst_durations);
    avg_ibi = sum(IBIs)/length(IBIs)
    std_ibi = std(IBIs)
    if verbose >= 1
        println("The average spike duration is $(round(avg_spike_dur, digits = 2)) ms +- $(round(std_spike_dur, digits = 2)) ms")
        println("The average burst duration is $(round(avg_burst_dur;digits = 2)) ms +- $(round(std_burst_dur;digits = 2)) ms")
        println("The average interburst interval is $(round(avg_ibi;digits = 2)) ms +- $(round(std_ibi;digits = 2)) ms")
    end
    if mode == 1
        return [avg_spike_dur, std_spike_dur, avg_burst_dur, std_burst_dur, avg_ibi, std_ibi]
    elseif mode == 2
        #Mode 2 returns lists of spike durations, bursts and ibis
        return spike_durations, burst_durations, IBIs
    end
end

"""
For 3D arrays and functions, this will extract all of the bursts and convert it into a graphable array
"""
function extract_burstmap(spike_array::BitArray{3})
    burst_map = zeros(size(spike_array))
    burst_data = max_interval_algorithim(spike_array; dt = 1.0)
    for (x,y,data) in burst_data
        if data[1] != nothing
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
