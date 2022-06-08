using Revise
using RetinalChaos

#%% Lets set up a SDE problem and then run it repeatedly for different values for gCa
conds_dict = read_JSON("params\\GABA_conds.json")
u0 = extract_dict(conds_dict, GABA_conds)#Initial conditions
pars_dict = read_JSON("params\\GABA_params.json")
pars_dict[:g_GABA] = 0.0
pars_dict[:g_ACh] = 0.0
p = extract_dict(pars_dict, GABA_pars) #Parameters
tspan = (0.0, 300e3) #Timespan
prob = SDEProblem(GABA_SDE, noise, u0, tspan, p) #ODE problem

#If we want to look at the baseline trace
@time sol = solve(prob, SOSRI());
plot(sol, vars=[1])

#%% Run the parameter space for gCa
n_trajectories = 10
par_idx = p_find(:g_Ca; list_p=GABA_pars) #Point to the index of the parameter
test_rng = LinRange(1.0, 20.0, n_trajectories) #Determine the range of the parameters (specified above)

prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #Set up the problem

@time sim = solve(ensemble_prob, SOSRI(), save_idxs=[1], trajectories=n_trajectories, EnsembleThreads(), progress=true, progress_steps=1);

# Plot some histograms
plt_a = plot()
plt_b1a = plot()
plt_b1b = plot()
plt_b2a = plot()
plt_b2b = plot()

for (sol_idx, sol_i) in enumerate(sim)
     println(test_rng[sol_idx])
     plot!(plt_a, sol_i, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
     #%% Plot each parameter in a histogram
     println(size(sol_i))
     thresh = calculate_threshold(sol_i)[1]
     spike_tstamps = get_timestamps(sol_i; threshold=thresh)
     spike_durs, isi = extract_interval(spike_tstamps)

     hfit = fit(Histogram, spike_durs, LinRange(5, 15, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b1a, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")

     hfit = fit(Histogram, isi, LinRange(1, 10, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b1b, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")

     burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
     burst_durs, ibi = extract_interval(burst_tstamps)

     hfit = fit(Histogram, burst_durs, LinRange(500, 2000, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b2a, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")

     hfit = fit(Histogram, ibi, LinRange(2000, 200000, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b2b, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")
end
plt_a
plt_b = plot(plt_b1a, plt_b1b, plt_b2a, plt_b2b, layout=(2, 2))
plt_CA = plot(plt_a, plt_b, layout=(2, 1))

#%% Run the parameter space for gK
n_trajectories = 10
par_idx = p_find(:g_K; list_p=GABA_pars) #Point to the index of the parameter
test_rng = LinRange(7.5, 12.5, n_trajectories) #Determine the range of the parameters (specified above)

prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #Set up the problem

@time sim = solve(ensemble_prob, SOSRI(), save_idxs =[1], trajectories=n_trajectories, EnsembleThreads(), progress=true, progress_steps=1);

# Plot some histograms
plt_a = plot()
plt_b1a = plot()
plt_b1b = plot()
plt_b2a = plot()
plt_b2b = plot()

for (sol_idx, sol_i) in enumerate(sim)
     println(test_rng[sol_idx])
     plot!(plt_a, sol_i, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
     #%% Plot each parameter in a histogram
     println(size(sol_i))
     thresh = calculate_threshold(sol_i)[1]
     spike_tstamps = get_timestamps(sol_i; threshold=thresh)
     spike_durs, isi = extract_interval(spike_tstamps)

     hfit = fit(Histogram, spike_durs, LinRange(1, 10, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b1a, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")

     hfit = fit(Histogram, isi, LinRange(1, 10, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b1b, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")

     burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
     burst_durs, ibi = extract_interval(burst_tstamps)

     hfit = fit(Histogram, burst_durs, LinRange(500, 2000, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b2a, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")

     hfit = fit(Histogram, ibi, LinRange(2000, 200000, 50))
     weights = hfit.weights / maximum(hfit.weights)
     edges = collect(hfit.edges[1])[1:length(hfit.weights)]
     plot!(plt_b2b, edges, weights, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), label="")
end
plt_a
plt_b = plot(plt_b1a, plt_b1b, plt_b2a, plt_b2b, layout=(2, 2))
plt_K = plot(plt_a, plt_b, layout=(2, 1))