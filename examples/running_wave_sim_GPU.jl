using Revise
using RetinalChaos
using Dates, Plots, BSON
using StatsBase, StatsPlots
#using CUDA, DifferentialEquations, ResettableStacks, RandomNumbers, LinearAlgebra
#Configure the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
import RetinalChaos.param_path

#Load the telegram client and env
dotenv("D:\\TelegramAccessEnv\\.env")
#Activate the GPU
RetinalChaos.CUDA.allowscalar(false)

#%% Run a longer simulation loop
#for mu in LinRange(0.05, 1.0, 25) #We want to rerun this exp with wave extraction
for mu in [0.0, 0.125, 0.18, 0.25, 0.50]
    println("Running simulation for $mu")
    try
        BotNotify("{Waves} Running simulation for 0% started")
        u_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/conds.json")
        p_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/params.json")
        p_dict[:μ] = 0.0
        save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\mu_experiment\\mu_$(round(Int64, mu*100))\\"
        NetSol = run_model(save_path, p_dict, u_dict)
        BotNotify("{Waves} Running simulation for 18% completed")
        save_solution(NetSol, save_path, partitions = 4)
        BotNotify("{Waves} Saving simulation for 18% completed")
        sol = load_solution(save_path)
        animate_solution(sol, save_path)
        BotNotify("{Waves} Animating simulation for 18% completed")
        timeseries_analysis(sol, save_path)
        BotNotify("{Waves} Timeseries analysis completed")
    catch error
        println(error)
        BotNotify("{Waves} An error has occurred $(typeof(error))")
    end
end

#%% Run a simulation on it's own
try
    BotNotify("{Waves} Running simulation for 18% started")
    u_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/conds.json") #Load initial conditions
    p_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/params.json") #Load the parameters
    ################# If you want to modify parameters do it here ###############
    p_dict[:t_run] = 300e3 #We want to make the run time longer
    p_dict[:μ] = 0.0
    ################# If you want to modify parameters do it here ###############
    save_path = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\example_animations\\mu_18" #Specify the save path
    NetSol = run_model(save_path, p_dict, u_dict) #This is a succinct function for running the model
    BotNotify("{Waves} Running simulation for 18% completed")
    save_solution(NetSol, save_path, partitions = 10) #Save the solution in 10 seperate parts
    BotNotify("{Waves} Saving simulation for 18% completed")
    sol = load_solution(save_path) #Reload the solution for analysis
    animate_solution(sol, save_path) #Animate the solution
    BotNotify("{Waves} Animating simulation for 18% completed")
    timeseries_analysis(sol, save_path) #conduct a timeseries analysis
    BotNotify("{Waves} Timeseries analysis completed")
catch 
    BotNotify("{Waves} An error has occurred $(typeof(error))")
end

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