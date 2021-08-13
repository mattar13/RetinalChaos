using Revise
using RetinalChaos #Load the package
using BSON, JLD2
#%%Lets figure out how to extract the waves
import RetinalChaos.param_path

load_path = "E:\\Data\\Modelling\\mu_12\\"
p_dict = read_JSON("$(load_path)\\params.json", is_type = Dict{Symbol, Float32})
u_dict = read_JSON("$(load_path)\\iconds.json", is_type = Dict{Symbol, Float32})
sol = load_model(load_path, p_dict, u_dict, gpu = false)
#%%

RetinalChaos.save_solution(sol, load_path)
#%%
RetinalChaos.load_solution(load_path)


#%%
timestamps = BSON.load("$(load_path)\\timestamps.bson") #Access timestamps
data = BSON.load("$(load_path)\\data.bson") #Access Thresholds
#%%
threshs = data["Thresholds"]
wave_markers = extract_waves(sol_mu12, threshs)
wave_tstamps = get_timestamps(wave_markers, sol_mu12.t)


spike_tstamps = timestamps["Spikes"]