#######################Fitting the models######################
using DiffEqBase



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

function IntervalFunc(Ŷ, p; verbose=true, mode = :All, kwargs...)

    conds_dict = read_JSON("params\\conds.json")
    u0 = conds_dict |> extract_dict
    tspan = (0.0, 300e3)
    prob = ODEProblem(T_ODE, u0, tspan, p)
    if verbose
        @time sol = solve(prob,  saveat=1.0) #So far the best method is SOSRI
    else
        sol = solve(prob, saveat=1.0) #So far the best method is SOSRI
    end
    timestamps = timeseries_analysis(sol.t, sol(sol.t, idxs=1), timestamps_only=true) #We can use a special version that only reports timestamps
    #println(timestamps)
    spike_durs, ISIs = extract_interval(timestamps["Spikes"])
    burst_durs, IBIs = extract_interval(timestamps["Bursts"])
    burst_dur = sum(burst_durs) / length(burst_durs)
    spike_dur = sum(spike_durs) / length(spike_durs)
    ISI = sum(ISIs) / length(ISIs)
    IBI = sum(IBIs) / length(IBIs)
    if mode == :SpikeDuration
        return spike_dur
    elseif mode == :SpikeInterval
        return ISI
    elseif mode == :BurstDuration 
        return burst_dur
    elseif mode == :BurstInterval
        return IBI
    elseif mode == :All
        return [spike_dur-Ŷ[1], ISI-Ŷ[2], burst_dur-Ŷ[3], IBI-Ŷ[4]]
    else
        throw("Incorrect setting")
    end
end