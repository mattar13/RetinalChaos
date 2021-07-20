using Revise
using RetinalChaos
using Dates, Plots, JLD2
using CUDA, DifferentialEquations, ResettableStacks, RandomNumbers, LinearAlgebra
using TimerOutputs #This is for benchmarking
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
u_dict = read_JSON("params\\conds.json", is_type = Dict{Symbol, Float32})

#%% Testing the conversion to CPU function

for mu in LinRange(0.05, 1.0, 25)
    #BotNotify("{Waves} Running simulation for mu = $mu")
    save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\mu_$(round(Int64, mu*100))\\"
    p_dict[:μ] = mu
    NetSol = load_model(save_path, p_dict, u_dict)
    timestamps, data = timeseries_analysis(save_path, NetSol)
    #Maybe we are running out of GPU memory and need to reset here
    NetSol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
    BotNotify("{Waves} Running simulation for mu = $mu completed")
end
#%% Load the BSON file for the timeseries_analysis

load_analysis = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\mu_5\\"
model = load_model(load_analysis, p_dict, u_dict)
#%%
threshold = calculate_threshold(model, Z = 6)
tstamps = get_timestamps(model, threshold)
spike_result = extract_interval(tstamps; max_duration = 100.0)
p1 = histogram(result[1], yaxis = :log)
p2 = histogram(result[2], yaxis = :log)
plot(p1, p2)
#%% Lets look at the bursts now

#%%
data = BSON.load("$(load_analysis)\\data.bson")





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