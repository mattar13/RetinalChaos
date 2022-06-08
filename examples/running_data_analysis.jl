# This file contains supplemental analysis figures
using Revise
using RetinalChaos
#Setup the fonts and stuff
import RetinalChaos: calculate_threshold, get_timestamps, max_interval_algorithim
#Step 1: Make and run the model
conds_dict = read_JSON("params\\conds.json")
conds_dict[:I_app] = 1.0
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
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps_var)
burst_durs, ibi = extract_interval(burst_tstamps)

#%% Lets plot things to see more
