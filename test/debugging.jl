using Revise
using RetinalChaos #Load the package
using BSON
using ImageSegmentation
#%%Lets figure out how to extract the waves
load_path = "E:\\Data\\Modelling\\mu_12\\"
p_dict = read_JSON("$(load_path)\\params.json", is_type = Dict{Symbol, Float32})
u_dict = read_JSON("$(load_path)\\iconds.json", is_type = Dict{Symbol, Float32})
sol_mu12 = load_model(load_path, p_dict, u_dict, gpu = false)
timestamps = BSON.load("$(load_path)\\timestamps.bson") #Access timestamps
data = BSON.load("$(load_path)\\data.bson") #Access Thresholds
#%%
threshs = data["Thresholds"]
extract_waves(sol_mu12, threshs)
#%%