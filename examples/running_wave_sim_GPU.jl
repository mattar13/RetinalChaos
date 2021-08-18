using Revise
using RetinalChaos
using Dates, Plots, JLD2
using StatsBase, StatsPlots
using CUDA, DifferentialEquations, ResettableStacks, RandomNumbers, LinearAlgebra
#using TimerOutputs #This is for benchmarking
#Configure the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())


#Load the telegram client and env
dotenv("D:\\TelegramAccessEnv\\.env")
#Activate the GPU
RetinalChaos.CUDA.allowscalar(false)

# Load the needed files to run the model
param_root = RetinalChaos.param_path
#params/params.json
p_dict = read_JSON(Dict{Symbol, Float32}, "params/params.json")
p_dict[:t_warm] = 1000.0
p_dict[:t_run] = 1000.0 #Extend the simulation time so we can find longer bursts
p_dict[:μ] = 0.18
u_dict = read_JSON(Dict{Symbol, Float32}, "params/conds.json")
save_path = "/home/john/Documents/modelling/mu_18/"
NetSol = run_model(save_path, p_dict, u_dict)
#%%
t = collect(1:1000) .* 5.0e-3

LinRange(1, length(t), 4)[1:end-1]


file_contents = Dict(
    :sol_t => ["sol_t.bson"],
    :sol_u => ["sol_u.bson"]
)
write_JSON(file_contents, "file_contents.json")
read_JSON(Dict{String, Array}, "file_contents.json")
#%%
#for mu in [0.0, 0.125, 0.25, 0.50]
for mu in LinRange(0.05, 1.0, 25) #We want to rerun this exp with wave extraction
    try
        BotNotify("{Waves} Running simulation for mu = $mu")
        save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\mu_$(round(Int64, mu*100))\\"
        p_dict[:μ] = mu
        NetSol = run_model(save_path, p_dict, u_dict)
        timestamps, data = timeseries_analysis(save_path, NetSol, reltol = 1e-2)
        #Maybe we are running out of GPU memory and need to reset here
        NetSol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
        println("{Waves} Running simulation for mu = $mu")
        BotNotify("{Waves} Running simulation for mu = $mu completed")
    catch error
        BotNotify("{Waves} has encountered an error: $error")
    end
end

#%% Run a simulation on it's own
#save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\mu_60"
#p_dict[:t_run] = 60e3 #Extend the simulation time so we can find longer bursts
#p_dict[:μ] = 0.6041667 #Change the parameter
#NetSol = load_model(save_path, p_dict, u_dict, abstol = 2e-2, reltol = 1e-2)
#timestamps, data = timeseries_analysis(save_path, NetSol)

#%% This is for a replication experiment
#for repeat in 1:4, mu in LinRange(0.05, 1.0, 25)
#    BotNotify("{Waves} Running simulation for mu = $mu")
#    save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\r_$(repeat)_mu_$(round(Int64, mu*100))\\"
#    p_dict[:μ] = mu
#    NetSol = load_model(save_path, p_dict, u_dict)
#    timestamps, data = timeseries_analysis(save_path, NetSol)
#    #Maybe we are running out of GPU memory and need to reset here
#    NetSol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
#    BotNotify("{Waves} Running simulation for mu = $mu completed")
#end