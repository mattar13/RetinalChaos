using Revise
using RetinalChaos, NeuroPhys

#%% Lets load a file from google drive
using GoogleDrive
drive_link = "https://drive.google.com/drive/folders/1AmUCVDjMbmg-zF_Me89F1PdV9Yxqt7C0?usp=sharing"
target_file = "E:\\Data\\Patching\\2019_11_03_Patch\\Animal_2\\Cell_3\\"
google_download(drive_link, target_file)



#%% Load some data
target_file = "E:\\Data\\Patching\\2019_11_03_Patch\\Animal_2\\Cell_3\\19n03042.abf"
data = extract_abf(target_file)
#Conduct the timescale analysis on this file
#reduce the size of the file
truncate_data!(data, 
     t_pre = 140.0, t_post = 200.0, 
     truncate_based_on = :time_range)
data-25
plot(data)
#%%

#%% Now simulate some traces
p = read_JSON("params/params.json") 
u0 = read_JSON("params/conds.json")
tspan = (0.0, 120e3);
prob_sde = SDEProblem(T_sde, noise, u0|>extract_dict, tspan, p|>extract_dict);
#Inject a current
#p[:I_app] = 0.42
prob_basic = ODEProblem(T_ode, u0|>extract_dict, tspan, p|>extract_dict)
sol_sde = solve(prob_sde, progress = true)
sol_basic = solve(prob_basic, progress = true)
#%%
plot(sol_sde, vars = [1])
plot!(sol_basic, vars = [1])
#%%
n = 10
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, :n, LinRange(0.0, 1.0, n))
ensemble_prob = EnsembleProblem(prob_sde, prob_func = prob_func);
@time sim = solve(ensemble_prob, trajectories = n, EnsembleThreads());

#%%
plt = plot()
for sol in sim
     plot!(plt, sol)
end
plt
