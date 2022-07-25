#######################Fitting the models######################
using DiffEqBase
"""
This function normalizes values to a range between values
"""
function normalize(a; old_rng = (nothing, nothing), rng = (0, 1))

    if isnothing(old_rng[1])
        min = minimum(a)
    else
        min = old_rng[1]
    end

    if isnothing(old_rng[2])
        max = maximum(a)
    else
        max = old_rng[2]
    end

    norm = max - min;
    normed = (a .- min) ./ norm;
    x, y = rng
    scale = y - x;
    (normed * scale) .+ x;
end

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
    sum/n
end


function MeanSquaredErrorSOL(Ŷ::T, sol; region=(-100, 2500)) where {T}
    timestamps = timeseries_analysis(sol.t, sol(t, idxs=1), timestamps_only=true) #We can use a special version that only reports timestamps
    println(timestamps)
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

function loss_func(Ŷ, p; region=(-100, 2500))
    conds_dict = read_JSON("params\\conds.json")
    u0 = conds_dict |> extract_dict
    tspan = (0.0, 300e3)
    prob = SDEProblem(T_SDE, noise, u0, tspan, p)
    @time sol = solve(prob, SOSRI(), saveat=1.0) #So far the best method is SOSRI
    MSE = MeanSquaredErrorSOL(Ŷ, sol)
    return MSE
    #vt = copy(sol(sol.t, idxs=1) |> Array)
    #timestamps = timeseries_analysis(sol.t, vt, timestamps_only = true) #We can use a special version that only reports timestamps
    #println(timestamps)
    #if !isnothing(timestamps["Spikes"])
    #    test_burst_idx = round.(Int64, timestamps["Spikes"][1, 1]+region[1]:1.0:timestamps["Spikes"][1, 1]+region[2])
    #    Y = sol(test_burst_idx, idxs = 1) |> Array
    #    MSE = MeanSquaredError(Ŷ, Y) #This is a point by point comparison between the traces
    #    return MSE
    #else
    #    return Inf #We really don't want these results
    #end
end