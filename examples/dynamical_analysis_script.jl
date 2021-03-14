#%% This shows all the different things you can do with the dynamical analysis of the model
using RetinalChaos
#Set up the plotting
font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)

#Set up the file root and default parameters
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

#save everything in the figures folder
save_figs = "figures\\"
if isdir(save_figs) == false
    #The directory does not exist, we have to make it 
    mkdir(save_figs)
end

#%% Supplemental figure, Bifurcation analysis of voltage injections
p = read_JSON(params_file) 
u0 = read_JSON(conds_file)
p[:I_app] = 0.0 #Set initial applied current to 0
p[:g_ACh] = 0.0 #Remove g_ACh influence
p[:g_TREK] = 0.0 #Remove g_TREK influence
tspan = (0.0, 30e3);
prob_eq = ODEProblem(T_ode, u0|>extract_dict, tspan, p|>extract_dict)

#%% Codim 1 analysis
codim1 = (:I_app)
c1_lims = (-40.0, 1.0)
print("Codimensional analysis time to complete:")
@time c1_map = codim_map(prob_eq, codim1, c1_lims = c1_lims, eq_res = 10)
plot(c1_map, xlabel = "Injected Current", ylabel = "Membrane Voltage")

#%% Codim 3 analysis
codim3 = (:g_K, :g_Ca, :I_app)
c1_lims = (1.0,20.0); c2_lims = (1.0, 20.0); c3_lims = (-50.0, 1.0) 
print("Codimensional analysis time to complete:")
@time c3_map = codim_map(prob_eq, codim3, c1_lims = c1_lims, c2_lims = c2_lims, resolution = 20);
c3_map.equilibria |> length
#%%
plot(c3_map, view = :zxy)

#%% Supplemental figure. Analysis of parameter change
test_rng = range(0.5, 2.0, length = 100) #this ranges from halving the parameter to doubling it
p = read_JSON(params_file); #Because of my catch, we can keep these as dictionaries 
u0 = read_JSON(conds_file);
#We can set the initial voltage sensitivty here 
p[:σ] = 10.0 
tspan = (0.0, 60e3)
prob = ODEProblem(T_ode, u0, tspan, p);
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
    println("Simulation for parameter $par completed")
    for (sol_idx, sol) in enumerate(sim)
        #println(sol |> length)
        dt = 0.1 #set the time differential according to supp figure 1
        t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
        v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
        ts_analysis = timescale_analysis(v_t, dt = dt)
        for i = 1:length(ts_analysis)
            results[par_idx, sol_idx, i] = sum(ts_analysis[i])/length(ts_analysis[i])
        end
    end
    println("Statistics for $par calculated")
end
#%%
for (idx, p) in enumerate(T_ode.ps)
    println("Plotting results for $p")
    sfig3i = plot(layout = grid(3,1), size = (1000, 800))
    plot!(sfig3i[1], test_rng, results[idx,:,1], label = "$p", xlabel = "Norm Par", ylabel = "Spike Duration")
    plot!(sfig3i[2], test_rng, results[idx,:,2]./1000, label = "$p", xlabel = "Norm Par", ylabel = "Burst Duration")
    plot!(sfig3i[3], test_rng, results[idx,:,3]./1000, label = "$p", xlabel = "Norm Par", ylabel = "IBI")
    savefig(sfig3i, joinpath(save_figs, "gradient_for_param_$p.png"))
end
#sfig3i