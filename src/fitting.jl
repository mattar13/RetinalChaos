#######################Fitting the models######################
"""
This calculated the Mean squared error or the L2 Norm
It takes in data the model has predicted (Ȳ) and experimentally observed data (Y)
usage:
    \n Ȳ = predicted data
    \n Y = observed data
    \n loss = MSE(Ŷ, Y)
"""
function MeanSquaredError(Ŷ, Y) where {T}
    n = length(Y)
    sum = 0.0
    for i in 1:n
        diff = Y[i] - Ŷ[i]
        squared_diff = diff^2
        sum += squared_diff
    end
    sum / n
end


function MeanSquaredErrorSOL(Ŷ::T, sol; region=(-100, 2500)) where {T}
    timestamps = timeseries_analysis(sol.t, sol(sol.t, idxs=1), timestamps_only=true) #We can use a special version that only reports timestamps
    #println(timestamps)
    if !isempty(timestamps["Spikes"])
        test_burst_idx = round.(Int64, timestamps["Spikes"][1, 1]+region[1]:1.0:timestamps["Spikes"][1, 1]+region[2])
        Y = sol(test_burst_idx, idxs=1)
        MSE = MeanSquaredError(Ŷ, Y) #This is a point by point comparison between the traces
        return MSE
    else
        return Inf
    end
end

"""
Can we write a function that creates a gradient? 
"""

function MSELoss(Ŷ, p; verbose=false, kwargs...)
    #println(p)
    conds_dict = read_JSON("params\\conds.json")
    u0 = conds_dict |> extract_dict
    tspan = (0.0, 300e3)
    prob = SDEProblem(T_SDE, noise, u0, tspan, p)
    if verbose
        @time sol = solve(prob, SOSRI(), saveat=1.0) #So far the best method is SOSRI
    else
        sol = solve(prob, SOSRI(), saveat=1.0) #So far the best method is SOSRI
    end
    MSE = MeanSquaredErrorSOL(Ŷ, sol; kwargs...)
    return MSE
end

"""
Can we write a function that creates a gradient? 
"""
#=
function extract_spike_trace(ts, data; dt = 1.0, cell_n = 1, idx=1, spike_dur=25, normalize = true)
    spike_begin = (ts["Spikes"][cell_n][idx, 1]/dt)
    #println(spike_begin)
    spike_rng = round.(Int64, spike_begin+1:spike_begin + Int64(spike_dur/dt))
    #print("Cell number $cell_n begins at: ")
    #println(spike_rng[1])
    t = data["Time"][spike_rng]

    spike = data["DataArray"][cell_n, spike_rng]
    if normalize
        spike = standardize(UnitRangeTransform, spike)
    end
    return t, spike
end

function extract_burst_trace(ts, data;  dt = 1.0, cell_n = 1, idx=1, burst_dur=1000, normalize = true)
    burst_begin = ts["Bursts"][cell_n][idx, 1]/dt
    burst_rng = round.(Int64, burst_begin+1:burst_begin+Int64(burst_dur/dt))
    t = data["Time"][burst_rng]
    burst = data["DataArray"][cell_n, burst_rng]
    if normalize
        burst = standardize(UnitRangeTransform, burst)
    end
    return t, burst
end

function extract_IBI_trace(ts, data; dt = 1.0, cell_n = 1, idx=1, IBI_dur=60e3, normalize = true)
    IBI_begin = ts["Bursts"][cell_n][idx, 1]/dt
    #println(IBI_begin)
    IBI_rng = round.(Int64, IBI_begin+1:IBI_begin + Int64(IBI_dur/dt))
    t = data["Time"][IBI_rng]
    IBI = data["DataArray"][cell_n, IBI_rng]
    if normalize
        IBI = standardize(UnitRangeTransform, IBI)
    end
    return t, IBI
end
=#

function extract_trace(ts::Dict{String,Matrix{T}}, data;
    tstamps="Spikes", dt=1.0, idx=1, duration=25, normalize=true
) where {T}
    trace_begin = (ts[tstamps][idx, 1] / dt)
    #println(spike_begin)
    trace_rng = round.(Int64, trace_begin+1:trace_begin+Int64(duration / dt))
    #print("Cell number $cell_n begins at: ")
    #println(spike_rng[1])
    t = data["Time"][trace_rng]

    trace = data["DataArray"][trace_rng]
    if normalize
        spike = standardize(UnitRangeTransform, trace)
    end
    return t, spike
end

function extract_trace(ts::Dict{String,Vector{Matrix{T}}}, data;
    tstamps="Spikes", dt=1.0, cell_n=1, idx=1, duration=25, normalize=true
) where {T}
    trace_begin = (ts[tstamps][cell_n][idx, 1] / dt)
    trace_rng = round.(Int64, trace_begin+1:trace_begin+Int64(duration / dt))
    t = data["Time"][trace_rng]
    trace = data["DataArray"][cell_n, trace_rng]
    if normalize
        trace = standardize(UnitRangeTransform, trace)
    end
    return t, trace
end

function extract_trace(ts::Dict{String,Matrix{Matrix{T}}}, data;
    tstamps="Spikes", dt=1.0, cell_swp=1, cell_ch = 1,  idx=1, duration=25, normalize=true
) where {T}
    trace_begin = (ts[tstamps][cell_swp, cell_ch][idx, 1] / dt)
    trace_rng = round.(Int64, trace_begin+1:trace_begin+Int64(duration / dt))
    t = data["Time"][trace_rng]
    trace = data["DataArray"][cell_swp, trace_rng, cell_ch]
    if normalize
        trace = standardize(UnitRangeTransform, trace)
    end
    return t, trace
end

extract_spike_trace(ts, data; spike_dur=25, kwargs...) = extract_trace(ts, data; tstamps="Spikes", duration=spike_dur, kwargs...)
extract_burst_trace(ts, data; burst_dur=1000, kwargs...) = extract_trace(ts, data; tstamps="Bursts", duration=burst_dur, kwargs...)
extract_IBI_trace(ts, data; IBI_dur=60e3, kwargs...) = extract_trace(ts, data; tstamps="Bursts", duration=IBI_dur, kwargs...)

#This is the case if only a sinle trace is being calculated against another single trace
function IntervalLoss(tsŶ::Dict{String,Matrix{T}}, dataŶ, tsY::Dict{String,Matrix{T}}, dataY;
    tstamps="Spikes", dt=1.0, duration=25, normalize=true
) where {T<:Real}
    MSE = 0.0
    len = 0
    nŶ = size(tsŶ[tstamps], 1)
    nY = size(tsY[tstamps], 1)
    for iy in 1:nY, iŷ in 1:nŶ
        #println(iy)
        #println(iŷ)
        if tsŶ[tstamps][iŷ, 1] + duration / dt > length(dataŶ["DataArray"])
            nothing
            #println(tsŶ[tstamps][iŷ, 1] + duration / dt)
        elseif tsY[tstamps][iy, 1] + duration / dt > length(dataY["DataArray"])
            nothing
            #println(tsY[tstamps][iy, 1] + duration / dt)
        else
            len += 1
            tŶ, spike_Ŷ = extract_trace(tsŶ, dataŶ, duration=duration, normalize=normalize, idx=iŷ)
            tY, spike_Y = extract_trace(tsY, dataY, duration=duration, normalize=normalize, idx=iy)
            MSE += MeanSquaredError(spike_Ŷ, spike_Y)
        end
        #println(len)
    end
    MSE / len
end

function IntervalLoss(tsŶ::Dict{String, Vector{Matrix{T}}}, dataŶ, tsY::Dict{String, Vector{Matrix{T}}}, dataY;
    dt=1.0, duration=25, tstamps="Spikes", normalize=true
) where {T<:Real}
    println("Vector of Matrix Version")
    MSE = 0.0
    len = 0
    #We have to do a different thing if the timestamps are from a single traces
    nCellsŶ = length(tsŶ[tstamps])
    nCellsY = length(tsY[tstamps])
    for nCellY in 1:nCellsY, nCellŶ in 1:nCellsŶ
        #println
        nŶ = size(tsŶ[tstamps][nCellŶ], 1)
        nY = size(tsY[tstamps][nCellY], 1)
        for iy in 1:nY, iŷ in 1:nŶ
            if tsŶ[tstamps][nCellŶ][iŷ, 1] + duration / dt > size(dataŶ["DataArray"], 2)
                nothing
                #println("Not enough space")
            elseif tsY[tstamps][nCellY][iy, 1] + duration / dt > size(dataY["DataArray"], 2)
                nothing
            else
                len += 1
                tŶ, spike_Ŷ = extract_trace(tsŶ, dataŶ, duration=duration, normalize=normalize, cell_n=nCellŶ, idx=iŷ)
                tY, spike_Y = extract_trace(tsY, dataY, duration=duration, normalize=normalize, cell_n=nCellY, idx=iy)
                MSE += MeanSquaredError(spike_Ŷ, spike_Y)
            end
        end
    end
    MSE / len
end

#This mode will compare an experiment to a 2D simulation
function IntervalLoss(tsŶ::Dict{String,Matrix{Matrix{T}}}, dataŶ, tsY::Dict{String,Vector{Matrix{T}}}, dataY;
    dt=1.0, duration=25, tstamps="Spikes", normalize=true
) where {T<:Real}
    #println("Vector of Matrix Version")
    MSE = 0.0
    len = 0
    #We have to do a different thing if the timestamps are from a single traces
    nCellsŶ_swp, nCellsŶ_ch = size(tsŶ[tstamps])
    nCellsY = length(tsY[tstamps])
    for nCellY in 1:nCellsY, nCellŶ_swp in 1:nCellsŶ_swp, nCellŶ_ch in 1:nCellsŶ_ch
        #println(nCellY, nCellŶ_swp, nCellŶ_ch)
        nŶ = size(tsŶ[tstamps][nCellŶ_swp, nCellŶ_ch], 1)
        nY = size(tsY[tstamps][nCellY], 1)
        for iy in 1:nY, iŷ in 1:nŶ
            if tsŶ[tstamps][nCellŶ_swp, nCellŶ_ch][iŷ, 1] + duration / dt > size(dataŶ["DataArray"], 2)
                nothing
                #println("Not enough space")
            elseif tsY[tstamps][nCellY][iy, 1] + duration / dt > size(dataY["DataArray"], 2)
                nothing
            else
                len += 1
                tŶ, spike_Ŷ = extract_trace(tsŶ, dataŶ, duration=duration, normalize=normalize, cell_swp=nCellŶ_swp, cell_ch = nCellŶ_ch, idx=iŷ)
                tY, spike_Y = extract_trace(tsY, dataY, duration=duration, normalize=normalize, cell_n=nCellY, idx=iy)
                MSE += MeanSquaredError(spike_Ŷ, spike_Y)
            end
        end
    end
    MSE / len
end

function IntervalLoss(tsŶ::Dict{String, Vector{Matrix{T}}}, dataŶ, tsY::Dict{String, Matrix{T}}, dataY;
    dt=1.0, duration=25, tstamps="Spikes", normalize=true
) where {T<:Real}
    #println("Vector of Matrix Version")
    MSE = 0.0
    len = 0
    #We have to do a different thing if the timestamps are from a single traces
    nCellsŶ = length(tsŶ[tstamps])
    for nCellŶ in 1:nCellsŶ
        #println
        nŶ = size(tsŶ[tstamps][nCellŶ], 1)
        nY = size(tsY[tstamps], 1)
        for iy in 1:nY, iŷ in 1:nŶ
            if tsŶ[tstamps][nCellŶ][iŷ, 1] + duration / dt > size(dataŶ["DataArray"], 2)
                nothing
                #println("Not enough space")
            elseif tsY[tstamps][iy, 1] + duration / dt > length(dataY["DataArray"])
                nothing
            else
                len += 1
                tŶ, spike_Ŷ = extract_trace(tsŶ, dataŶ, duration=duration, normalize=normalize, cell_n=nCellŶ, idx=iŷ)
                tY, spike_Y = extract_trace(tsY, dataY, duration=duration, normalize=normalize, idx=iy)
                MSE += MeanSquaredError(spike_Ŷ, spike_Y)
            end
        end
    end
    MSE / len
end

#This mode will compare an experiment to a 1D simulation
function IntervalLoss(tsŶ::Dict{String, Matrix{Matrix{T}}}, dataŶ, tsY::Dict{String, Matrix{T}}, dataY;
    dt=1.0, duration=25, tstamps="Spikes", normalize=true
) where {T<:Real}
    MSE = 0.0
    len = 0
    #We have to do a different thing if the timestamps are from a single traces
    nCellsŶ_swp, nCellsŶ_ch = size(tsŶ[tstamps])
    for nCellŶ_swp in 1:nCellsŶ_swp, nCellŶ_ch in 1:nCellsŶ_ch
        #println
        nŶ = size(tsŶ[tstamps][nCellŶ_swp, nCellŶ_ch], 1)
        nY = size(tsY[tstamps], 1)
        for iy in 1:nY, iŷ in 1:nŶ
            if tsŶ[tstamps][nCellŶ_swp, nCellŶ_ch][iŷ, 1] + duration / dt > size(dataŶ["DataArray"], 2)
                nothing
                #println("Not enough space")
            elseif tsY[tstamps][iy, 1] + duration / dt > length(dataY["DataArray"])
                nothing
            else
                len += 1
                tŶ, spike_Ŷ = extract_trace(tsŶ, dataŶ, duration=duration, normalize=normalize, cell_swp=nCellŶ_swp, cell_ch=nCellŶ_ch, idx=iŷ)
                tY, spike_Y = extract_trace(tsY, dataY, duration=duration, normalize=normalize, idx=iy)
                MSE += MeanSquaredError(spike_Ŷ, spike_Y)
            end
        end
    end
    MSE / len
end

function IntervalLoss(tsŶ::Dict{String, Matrix{T}}, dataŶ, tsY::Dict{String, Vector{Matrix{T}}}, dataY;
    dt=1.0, duration=25, tstamps="Spikes", normalize=true
) where {T<:Real}
    #println("Vector of Matrix Version")
    MSE = 0.0
    len = 0
    #We have to do a different thing if the timestamps are from a single traces
    nCellsY = length(tsY[tstamps])
    for nCellY in 1:nCellsY
        #println
        nŶ = size(tsŶ[tstamps], 1)
        nY = size(tsY[tstamps][nCellY], 1)
        for iy in 1:nY, iŷ in 1:nŶ
            if tsŶ[tstamps][nCellŶ][iŷ, 1] + duration / dt > length(dataŶ["DataArray"])
                nothing
                #println("Not enough space")
            elseif tsY[tstamps][nCellY][iy, 1] + duration / dt > size(dataY["DataArray"], 2)
                nothing
            else
                len += 1
                tŶ, spike_Ŷ = extract_trace(tsŶ, dataŶ, duration=duration, normalize=normalize, idx=iŷ)
                tY, spike_Y = extract_trace(tsY, dataY, duration=duration, normalize=normalize, cell_n=nCellY, idx=iy)
                MSE += MeanSquaredError(spike_Ŷ, spike_Y)
            end
        end
    end
    MSE / len
end

function TimescaleLoss(tsŶ, dataŶ, tsY, dataY;
    dt=1.0, spike_dur=25, burst_dur=1000, IBI_dur=60e3,
    spike_weight=1.0, burst_weight=1.0, IBI_weight=1.0,
    normalize=true
)
    #We can change the data to normalize
    spike_MSE = IntervalLoss(tsŶ, dataŶ, tsY, dataY; dt=dt, tstamps = "Spikes", duration=spike_dur, normalize=normalize) * spike_weight
    #println(spike_MSE)
    burst_MSE = IntervalLoss(tsŶ, dataŶ, tsY, dataY; dt=dt, tstamps = "Bursts", duration=burst_dur, normalize=normalize) * burst_weight
    #println(burst_MSE)
    IBI_MSE = IntervalLoss(tsŶ, dataŶ, tsY, dataY; dt=dt, tstamps = "Bursts", duration=IBI_dur, normalize=normalize) * IBI_weight
    #println(IBI_MSE)
    #We can either return each number individually or sum them
    return spike_MSE + burst_MSE + IBI_MSE
end

function TimescaleLoss(Ŷ, Y; dt=1.0, kwargs...)
    tsŶ, dataŶ = timeseries_analysis(Ŷ)
    tsY, dataY = timeseries_analysis(Y)
    TimescaleLoss(tsŶ, dataŶ, tsY, dataY; dt=dt, kwargs...)
end

#This will allow us to calculate loss on a single parameter set
function TimescaleLoss(Ŷ, p::Vector{T}; idxs = 1, kwargs...) where {T<:Real}
    dt = Ŷ.dt
    conds_dict = read_JSON("params\\conds.json")
    u0 = conds_dict |> extract_dict
    tspan = (0.0, 300e3)
    prob = SDEProblem(T_SDE, noise, u0, tspan, p)
    sol = solve(prob, SOSRI(), save_idxs = idxs, saveat=dt, progress=true, progress_steps=1)
    TimescaleLoss(Ŷ, sol; dt=dt, kwargs...)
end
