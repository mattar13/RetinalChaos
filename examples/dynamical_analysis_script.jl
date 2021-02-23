using RetinalChaos
#%% This shows all the different things you can do with the dynamical analysis of the model

#%% Supplemental figure, Bifurcation analysis of voltage injections
p_dict = read_JSON(params_file) 
p_dict[:I_app] = 0.0 #Set initial applied current to 0
p_dict[:g_ACh] = 0.0 #Remove g_ACh influence
p_dict[:g_TREK] = 0.0 #Remove g_TREK influence
p = p_dict |> extract_dict;
u0 = read_JSON(conds_file) |> extract_dict;
tspan = (0.0, 30e3);
prob_eq = ODEProblem(T_ode, u0, tspan, p)
codim1 = (:I_app)
c1_lims = (-50.0, 20.0)
c1_map = codim_map(prob_eq, codim1, c1_lims = c1_lims, eq_res = 10)
plot(c1_map, xlabel = "Injected Current", ylabel = "Membrane Voltage")

#%% Supplemental figure. Analysis of parameter change
test_rng = range(0.5, 2.0, length = 50) #this ranges from halving the parameter to doubling it
p_dict = read_JSON(params_file) 
p_dict[:I_app] = 10.0
#We can set the initial voltage sensitivty here 
p = p_dict |> extract_dict;
u0 = read_JSON(conds_file) |> extract_dict;
tspan = (0.0, 60e3)
prob = ODEProblem(T_ode, u0, tspan, p_dict |> extract_dict);
results = zeros(length(T_ode.ps), length(test_rng), 3)
#Walk through each parameter
for (par_idx, par) in enumerate(T_ode.ps)
    println("Simulating test range for parameters: $par")
    #Setup an ensemble function to test each parameter
    par_rng = test_rng * p_dict[Symbol(par)]
    prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, par_rng)
    ensemble_prob = EnsembleProblem(prob, prob_func = prob_func);
    
    print("Time it took to simulate $(tspan[2]/1000)s:")
    @time sim = solve(ensemble_prob, EnsembleThreads(), trajectories = length(test_rng)); 
    for (sol_idx, sol) in enumerate(sim)
        println(sol |> length)
        dt = 0.1 #set the time differential according to supp figure 1
        t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
        v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
        ts_analysis = timescale_analysis(v_t, dt = dt)
        for i = 1:length(ts_analysis)
            results[par_idx, sol_idx, i] = sum(ts_analysis[i])/length(ts_analysis[i])
        end
    end
end
#%%
for (idx, p) in enumerate(T_ode.ps)
    println("Plotting results for $p")
    sfig3i = plot(layout = grid(3,1), size = (1000, 800))
    plot!(sfig3i[1], test_rng.*p_dict[Symbol(p)], results[idx,:,1], label = "$p")
    plot!(sfig3i[2], test_rng.*p_dict[Symbol(p)], results[idx,:,2]./1000, label = "$p")
    plot!(sfig3i[3], test_rng.*p_dict[Symbol(p)], results[idx,:,3]./1000, label = "$p")
    savefig(sfig3i, joinpath(save_figs, "gradient_for_param_$p.png"))
end
sfig3i