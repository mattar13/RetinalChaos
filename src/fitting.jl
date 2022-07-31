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
function extract_spike_trace(ts, data; dt = 1.0, cell_n = 1, idx=1, spike_dur=25, normalize = true)
    spike_begin = ts["Spikes"][cell_n][idx, 1]/dt
    spike_rng = round.(Int64, spike_begin:spike_begin + Int64(spike_dur/dt))
    t = data["Time"][spike_rng]
    spike = data["DataArray"][cell_n, spike_rng]
    if normalize
        spike = standardize(UnitRangeTransform, spike)
    end
    return t, spike
end

function extract_burst_trace(ts, data;  dt = 1.0, cell_n = 1, idx=1, burst_dur=1000, normalize = true)
    burst_begin = ts["Bursts"][cell_n][idx, 1]/dt
    burst_rng = round.(Int64, burst_begin:burst_begin+Int64(burst_dur/dt))
    t = data["Time"][burst_rng]
    burst = data["DataArray"][cell_n, burst_rng]
    if normalize
        burst = standardize(UnitRangeTransform, burst)
    end
    return t, burst
end

function extract_IBI_trace(ts, data; dt = 1.0, cell_n = 1, idx=1, IBI_dur=60e3, normalize = true)
    IBI_begin = ts["Bursts"][cell_n][idx, 1]/dt
    println(IBI_begin)
    IBI_rng = round.(Int64, IBI_begin:IBI_begin + Int64(IBI_dur/dt))
    t = data["Time"][IBI_rng]
    IBI = data["DataArray"][cell_n, IBI_rng]
    if normalize
        IBI = standardize(UnitRangeTransform, IBI)
    end
    return t, IBI
end

function SpikeLoss(tsŶ, dataŶ, tsY, dataY; dt=1.0, spike_dur=25, normalize = true)
    nŶ = size(tsŶ["Spikes"][1], 1)
    nY = size(tsY["Spikes"][1], 1)
    MSE = 0.0
    len = 0
    for iy in 1:nY, iŷ in 1:nŶ
        if tsŶ["Spikes"][1][iŷ, 1] + spike_dur/dt > size(dataŶ["DataArray"], 2)
            nothing
            #println("Not enough space")
        elseif tsY["Spikes"][1][iy, 1] + spike_dur/dt > size(dataY["DataArray"], 2)
            nothing
        else
            len += 1
            #println("Also not enough length")
            tŶ, spike_Ŷ = extract_spike_trace(tsŶ, dataŶ, idx=iŷ, spike_dur=spike_dur, normalize = normalize)
            tY, spike_Y = extract_spike_trace(tsY, dataY, idx=iy, spike_dur=spike_dur, normalize = normalize)
            MSE += MeanSquaredError(spike_Ŷ, spike_Y)
        end
    end    
    MSE / len
end

function BurstLoss(tsŶ, dataŶ, tsY, dataY; dt=1.0, burst_dur=1000, normalize = true)
    nŶ = size(tsŶ["Bursts"][1], 1)
    nY = size(tsY["Bursts"][1], 1)
    MSE = 0.0
    len = 0
    for iy in 1:nY, iŷ in 1:nŶ
        if tsŶ["Bursts"][1][iŷ, 1] + burst_dur/dt > size(dataŶ["DataArray"], 2)
            nothing
            #println("Not enough space")
        elseif tsY["Bursts"][1][iy, 1] + burst_dur/dt > size(dataY["DataArray"], 2)
            nothing
            #println("Also not enough length")
        else
            len += 1
            tŶ, burst_Ŷ = extract_burst_trace(tsŶ, dataŶ; idx=iŷ, burst_dur=burst_dur, normalize = normalize)
            tY, burst_Y = extract_burst_trace(tsY, dataY; idx=iy, burst_dur=burst_dur, normalize = normalize)
            MSE += MeanSquaredError(burst_Ŷ, burst_Y)
        end
    end
    MSE/len
end

function IBILoss(tsŶ, dataŶ, tsY, dataY; dt=1.0, IBI_dur = 60e3, normalize = true)
    nŶ = size(tsŶ["Bursts"][1], 1)
    nY = size(tsY["Bursts"][1], 1)
    MSE = 0.0
    len = 0
    for iy in 1:nY, iŷ in 1:nŶ
        if tsŶ["Bursts"][1][iŷ, 1] + IBI_dur/dt > size(dataŶ["DataArray"], 2)
            nothing
            #println("Not enough space")
        elseif tsY["Bursts"][1][iy, 1] + IBI_dur/dt > size(dataY["DataArray"], 2)
            nothing
            #println("Also not enough length")
        else
            len += 1
            tŶ, IBI_Ŷ = extract_IBI_trace(tsŶ, dataŶ; idx=iŷ, IBI_dur=IBI_dur, normalize = normalize)
            tY, IBI_Y = extract_IBI_trace(tsY, dataY; idx=iy, IBI_dur=IBI_dur, normalize = normalize)
            MSE += MeanSquaredError(IBI_Ŷ, IBI_Y)
        end
    end
    MSE / len
end

function TimescaleLoss(tsŶ, dataŶ, tsY, dataY;
    dt=1.0, spike_dur=25, burst_dur=1000, IBI_dur=60e3,
    spike_weight=1.0, burst_weight=1.0, IBI_weight=1.0, 
    normalize = true
)
    #We can change the data to normalize
    spike_MSE = SpikeLoss(tsŶ, dataŶ, tsY, dataY; dt = dt, spike_dur = spike_dur, normalize = normalize) * spike_weight
    println(spike_MSE)
    burst_MSE = BurstLoss(tsŶ, dataŶ, tsY, dataY; dt=dt, burst_dur = burst_dur, normalize = normalize) * burst_weight
    println(burst_MSE)
    IBI_MSE = IBILoss(tsŶ, dataŶ, tsY, dataY; dt=dt, IBI_dur=IBI_dur, normalize = normalize) * IBI_weight
    println(IBI_MSE)
    #We can either return each number individually or sum them
    return spike_MSE + burst_MSE + IBI_MSE
end

function TimescaleLoss(Ŷ, Y; dt = 1.0, kwargs...)
    tsŶ, dataŶ = timeseries_analysis(Ŷ)
    tsY, dataY = timeseries_analysis(Y)
    TimescaleLoss(tsŶ, dataŶ, tsY, dataY; dt = dt, kwargs...)
end

#This will allow us to calculate loss on a single parameter set
function TimescaleLoss(Ŷ, p::Vector{T}; kwargs...) where T <: Real
    dt = Ŷ.dt
    conds_dict = read_JSON("params\\conds.json")
    u0 = conds_dict |> extract_dict
    tspan = (0.0, 300e3)
    prob = SDEProblem(T_SDE, noise, u0, tspan, p)
    sol = solve(prob, SOSRI(), saveat=dt, progress=true, progress_steps=1); 
    TimescaleLoss(Ŷ, sol; dt = dt, kwargs...)
end
