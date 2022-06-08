using Revise
using RetinalChaos

#%% Lets set up a SDE problem and then run it repeatedly for different values for gCa
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict #Initial conditions
pars_dict = read_JSON("params\\params.json")
p = pars_dict |> extract_dict #Parameters
tspan = (0.0, 300e3) #Timespan
prob = SDEProblem(T_ODE, noise, u0, tspan, p) #ODE problem

n_trajectories = 10
par_idx = p_find(:g_Ca; list_p=GABA_pars) #Point to the index of the parameter
test_rng = LinRange(1.0, 20.0, n_trajectories) #Determine the range of the parameters (specified above)

prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #Set up the problem

@time sim = solve(ensemble_prob, trajectories=n_trajectories, EnsembleThreads());

plt_a = plot(sim[1], vars=[1], c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
for (sol_idx, sol_i) in enumerate(sim)
     println(test_rng[sol_idx])
     plot!(plt_a, sol_i, vars=[1], c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
end
plt_a