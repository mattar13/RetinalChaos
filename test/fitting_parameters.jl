using Revise
using RetinalChaos, NeuroPhys

#These parameters are the most likely to need to be fitted: 
#1) gGABA = [0.7 -> 1.0]
#2) gACh = [0.1 -> 0.215]
#2) gCa = [7.5 -> 8.0]
#3) gK = [3.0 -> 5.0]
#4) μ = [0.7-0.8]

#How does each parameter change the dataspace


#%% Supplemental figure. Analysis of parameter change
test_rng = range(0.5, 2.0, length=100) #this ranges from halving the parameter to doubling it
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
     ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

     print("Time it took to simulate $(tspan[2]/1000)s:")
     @time sim = solve(ensemble_prob, EnsembleThreads(), trajectories=length(test_rng))
     println("Simulation for parameter $par completed")
     for (sol_idx, sol) in enumerate(sim)
          #println(sol |> length)
          dt = 0.1 #set the time differential according to supp figure 1
          t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
          v_t = map(t -> sol(t)[1], t_rng) #extract according to the interval
          ts_analysis = timescale_analysis(v_t, dt=dt)
          for i = 1:length(ts_analysis)
               results[par_idx, sol_idx, i] = sum(ts_analysis[i]) / length(ts_analysis[i])
          end
     end
     println("Statistics for $par calculated")
end
#%%
for (idx, p) in enumerate(T_ode.ps)
     println("Plotting results for $p")
     sfig3i = plot(layout=grid(3, 1), size=(1000, 800))
     plot!(sfig3i[1], test_rng, results[idx, :, 1], label="$p", xlabel="Norm Par", ylabel="Spike Duration")
     plot!(sfig3i[2], test_rng, results[idx, :, 2] ./ 1000, label="$p", xlabel="Norm Par", ylabel="Burst Duration")
     plot!(sfig3i[3], test_rng, results[idx, :, 3] ./ 1000, label="$p", xlabel="Norm Par", ylabel="IBI")
     savefig(sfig3i, joinpath(save_figs, "gradient_for_param_$p.png"))
end
#sfig3i

