using Revise
using RetinalChaos #Load the package
using BSON
#%%Lets figure out how to extract the waves
import RetinalChaos.param_path

load_path = "E:\\Data\\Modelling\\mu_12\\"
p_dict = read_JSON("$(load_path)\\params.json", is_type = Dict{Symbol, Float32})
u_dict = read_JSON("$(load_path)\\iconds.json", is_type = Dict{Symbol, Float32})
sol_mu12 = load_model(load_path, p_dict, u_dict, gpu = false)
#%%
sol = sol_mu12
new_sol = SciMLBase.build_solution(sol.prob, sol.alg, sol.t, sol.u) #we can use this to build a solution without GPU
#%%
#Save everything but the solution
bson("$(load_path)\\testing_saving.bson", 
          Dict(:sol_prob => sol.prob, :sol_alg => sol.alg, :sol_t => sol.t)
     )

#%%
timestamps = BSON.load("$(load_path)\\timestamps.bson") #Access timestamps
data = BSON.load("$(load_path)\\data.bson") #Access Thresholds
#%%
threshs = data["Thresholds"]
wave_markers = extract_waves(sol_mu12, threshs)
wave_tstamps = get_timestamps(wave_markers, sol_mu12.t)


spike_tstamps = timestamps["Spikes"]