using Revise
using RetinalChaos #Load the package
using BSON
#%%Lets figure out how to extract the waves
import RetinalChaos.param_path

load_path = "E:\\Data\\Modelling\\mu_12\\"
p_dict = read_JSON("$(load_path)\\params.json", is_type = Dict{Symbol, Float32})
u_dict = read_JSON("$(load_path)\\iconds.json", is_type = Dict{Symbol, Float32})
sol = load_model(load_path, p_dict, u_dict, gpu = false)
#%%

RetinalChaos.save_solution(sol, load_path)
#%%
sol_data = BSON.load("$(load_path)\\sol_data.bson")
#%%
#remake net
net = Network()
#%%
bson("$(load_path)\\sol_2.bson", 
     Dict(:sol_d => Array(sol))
)
#%% Loading solution
reconstruct_u = map(x -> Vector(sol[:, x]), 1:length(sol.t))
sol_load = SciMLBase.build_solution(sol.prob, sol.alg, sol.t, reconstruct_u) #we can use this to build a solution without GPU
using StochasticDiffEq, LinearAlgebra

#%%
timestamps = BSON.load("$(load_path)\\timestamps.bson") #Access timestamps
data = BSON.load("$(load_path)\\data.bson") #Access Thresholds
#%%
threshs = data["Thresholds"]
wave_markers = extract_waves(sol_mu12, threshs)
wave_tstamps = get_timestamps(wave_markers, sol_mu12.t)


spike_tstamps = timestamps["Spikes"]