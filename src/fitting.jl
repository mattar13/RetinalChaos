#######################Fitting the models######################

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
function MSE(Ŷ, Y)
    n = length(Y)
    sum = 0.0
    for i in 1:n
        diff = Y[i] - Ŷ[i]
        squared_diff = diff^2
        sum += squared_diff
    end
    sum/n
end

"""
Data should be normalized before using this function
The accuracy is to the nearest x milliseconds
"""
function lagged_loss(Ŷ, Y; accuracy = 1000, weights = [1, 1])
    t_measure = collect(0:(length(Ŷ)-length(Y))/accuracy)
    loss_list = zeros(length(t_measure))
    for (idx, t0) in enumerate(t_measure)
        t_idx = t0 == 0 ? 1 : t0*1000
        fit_rng = Int.(collect(t_idx:(t_idx+length(Y)-1)))
        if any(Ŷ[fit_rng] .> 0.5)
            loss = MSE(Ŷ[fit_rng], Y)
            line_loss = sum(abs.((Ŷ[fit_rng] .> 0.5) .- (Y .> 0.5)))/1000
            loss_list[idx] = (weights[1]*loss) + (weights[2]*line_loss)
        else
            loss_list[idx] = Inf
        end
    end
    lag = t_measure[argmin(loss_list)] == 0.0 ? 1 : t_measure[argmin(loss_list)]  * 1000
    minimum(loss_list), lag
end


function monte_func(prob, i, repeat; pars = :all, dists = [Normal()])
    if pars == :all
        new_pars = map(x -> typeof(x) != Float64 ? rand(x) : x , dists);
        if new_pars[(:σ |> p_find)[1]] > 0.0
            stochastic = true
        end
        #new_conds = map(x -> typeof(x) != Float64 ? rand(x) : x,  dists);
        if stochastic
            return SDEProblem(BurstModel, noise, prob.u0, prob.tspan, new_pars)
        else
            return ODEProblem(BurstModel, prob.u0, prob.tspan, new_pars)
        end
    else
        stochastic = false
        new_pars = prob.p
        new_conds = prob.u0
        if prob.p[(:σ |> p_find)[1]] > 0.0
            stochastic = true
        end
        for (idx, var) in enumerate(pars)
            conds = var |> u_find
            pars = var |> p_find
            if length(conds) > 0
                new_conds[conds[1]] = rand(dists[idx])
            elseif length(pars) > 0
                new_pars[pars[1]] = rand(dists[idx])
            end
        end

        if stochastic
            return SDEProblem(BurstModel, noise, new_conds, prob.tspan, new_pars)
        else
            return ODEProblem(BurstModel, new_conds, prob.tspan, new_pars)
        end
    end
end

function evolve(prob, d_dict, Y_norm; 
        alg = SOSRI(), abstol = 2e-1, reltol = 2e-1, maxiters = 1e7, saveat = 1.0,
        loss_fxn = lagged_loss(),  #Arguments for weights for loss calculation
        par = :all, n_sims = 5, δ = 10, 
        iterations = 5, error = 5,
        verbose = false #Arguments for running the Monte Carlo simulation
    )
    loss_record = []
    old_rng = (minimum(Y), maximum(Y))
    dist_record = copy(d_dict)
    #Calculate first loss
    optimal_sim = prob
    optimal_p = prob.p
    dParLoss = zeros(length(prob.p), iterations)
    it = 0
    if verbose
        println("Simulating initial point:")
    end
    
    @time sol = solve(optimal_sim, abstol = abstol, reltol = reltol, maxiters = maxiters, saveat = saveat);
    Ŷ = map(t -> sol(t)[1], sol.t);Ŷ_norm = Ŷ |> f_norm; #Normalize initial values
    min_loss, lag = loss_fxn(Ŷ_norm, Y_norm) #calculate initial loss
    if verbose
        println("Initial Loss: $(min_loss)")
    end
    while it < iterations
        it += 1
        d0 = extract_dict(dist_record, BurstModel.params);
        pf = (pr, i, repeat) -> monte_func(pr, i, repeat; pars = par, dists = d0)
        ensemble_prob = EnsembleProblem(optimal_sim, prob_func = pf)
        if verbose
            println("Iteration $(it): Running Parallell simulations:")
        end
        @time sim = solve(ensemble_prob, alg, 
            abstol = abstol, reltol = reltol, maxiters = maxiters, 
            saveat = saveat, trajectories = n_sims, EnsembleThreads());
        if verbose
            println("Iteration $(it): Calculating losses")
        end
        losses = []
        for traj in sim
            Ŷ = map(t -> traj(t)[1], traj.t)
            Ŷ_norm = Ŷ |> x -> normalize(x; old_rng = old_rng)
            try
                loss, _ = lagged_loss(Ŷ_norm, Y_norm)
                push!(losses, loss)
            catch
                println("A parameter caused either a instability or maximum iterations")
                push!(losses, Inf)
            end
        end
        #Find out the minimal loss. If this is less that then initial loss
        if minimum(losses) < min_loss
            
            if verbose
                println("Iteration $(it): Loss has improved $(min_loss) -> $(minimum(losses))")
            end
            
            min_loss = minimum(losses)
            optimal_sim = sim[argmin(losses)].prob
            optimal_p = optimal_sim.p
            dParLoss[:,it] = optimal_p
            dist_record = create_distributions(BurstModel.params, optimal_p)
        else
            error -= 1
            dParLoss[:,it] = optimal_p
            n_sims += δ
            if verbose
                println("Iteration $(it): Loss has not improved $(min_loss) <=> $(minimum(losses))")
                println("Iteration $(it+1): Will run $(n_sims)")
            end
        end
        push!(loss_record, min_loss)
    end
    return optimal_sim, loss_record, dParLoss
end