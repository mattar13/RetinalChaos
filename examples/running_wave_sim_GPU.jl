using Revise
using RetinalChaos
using Dates, Plots, JLD2
using CUDA, DifferentialEquations, ResettableStacks, RandomNumbers, LinearAlgebra
#Configure the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#Load the telegram client and env
dotenv("D:\\TelegramAccessEnv\\.env")
#Activate the GPU
RetinalChaos.CUDA.allowscalar(false)

# Load the needed files to run the model
p_dict = read_JSON("params\\params.json", is_type = Dict{Symbol, Float32})
p_dict[:t_warm] = 10.0
p_dict[:t_run] = 10.0
u_dict = read_JSON("params\\conds.json", is_type = Dict{Symbol, Float32})

#%%
load_test = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\test_1"
model = load_model(load_test, p_dict, u_dict, animate_solution = false)



#%%
using BSON
#%%
BSON.@load "$(load_test)\\conds.bson" warmup #This is the JSON method
#%%
BSON.@load "$(load_test)\\sol.bson" sol
#%%
model = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
#%% 
#for mu in LinRange(0.05, 1.0, 25)
for mu in [0.88125,0.920833,0.960417,1.0]
    BotNotify("{Waves} Running simulation for mu = $mu")
    save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\mu_$(round(Int64, mu*100))\\"
    p_dict[:μ] = mu
    NetSol = load_model(save_path, p_dict, u_dict)
    timestamps, data = timeseries_analysis(save_path, NetSol)
    #Maybe we are running out of GPU memory and need to reset here
    NetSol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
    #BotNotify("{Waves} Running simulation for mu = $mu completed")
end

#%%
for repeat in 1:4, mu in LinRange(0.05, 1.0, 25)
    BotNotify("{Waves} Running simulation for mu = $mu")
    save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\r_$(repeat)_mu_$(round(Int64, mu*100))\\"
    p_dict[:μ] = mu
    NetSol = load_model(save_path, p_dict, u_dict)
    timestamps, data = timeseries_analysis(save_path, NetSol)
    #Maybe we are running out of GPU memory and need to reset here
    NetSol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
    BotNotify("{Waves} Running simulation for mu = $mu completed")
end