# This file contains supplemental analysis figures
using Revise
using RetinalChaos
#Setup the fonts and stuff
import RetinalChaos: calculate_threshold, get_timestamps, max_interval_algorithim
#Step 1: Make and run the model
conds_dict = read_JSON("params\\conds.json")
#conds_dict[:I_app] = 1.0
u0 = conds_dict |> extract_dict #Conditions
pars_dict = read_JSON("params\\params.json")
p = pars_dict |> extract_dict #Parameters
tspan = (0.0, 300e3) #Timespan
prob = SDEProblem(T_SDE, noise, u0, tspan, p) #Make the problem
@time sol = solve(prob, save_idxs=[1]); #Solve
plot(sol)

#%% Step 2: Conduct the analysis
thresh = calculate_threshold(sol)[1]
spike_tstamps = get_timestamps(sol; threshold=thresh)
spike_durs, isi = extract_interval(spike_tstamps)
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
burst_durs, ibi = extract_interval(burst_tstamps)

#%% Lets plot things to see more
#Zoom in on the first burst
xlims = (burst_tstamps[1][2, 1] - 50, burst_tstamps[1][2, 1] + 200)
plt_a = plot(sol)
vline!(plt_a, spike_tstamps[1][:, 1], c=:green, label="spike begin")
vline!(plt_a, spike_tstamps[1][:, 2], c=:red, label="spike end", xlims=xlims)

#Zoom in on the Bursts
plt_b = plot(sol)
vline!(plt_b, burst_tstamps[1][:, 1], c=:green, label="burst begin")
vline!(plt_b, burst_tstamps[1][:, 2], c=:red, label="burst end")

plt_b