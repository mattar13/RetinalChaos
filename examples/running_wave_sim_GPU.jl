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
p_dict = read_JSON("params/params.json", is_type = Dict{Symbol, Float32})
p_dict[:t_run] = 100e3 #Extend the simulation time so we can find longer bursts
p_dict[:μ] = 0.18
u_dict = read_JSON("params/conds.json", is_type = Dict{Symbol, Float32})
save_path = "/home/john/Documents/modelling/mu_18"

u0 = extract_dict(u_dict, p_dict[:nx], p_dict[:ny]) |> cu
p0 = p_dict |> extract_dict

#NetSol = run_model(save_path, p_dict, u_dict)
net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ], gpu = true)
NetProb = SDEProblem(net, noise, u0, (0f0 , p_dict[:t_warm]), p0)

print("[$(now())]: Warming up the solution... ")
@time sol = solve(NetProb, SOSRI(), 
    abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
    progress = true, progress_steps = 1, 
    save_everystep = false
)

warmup = sol[end]

NetProb = SDEProblem(net, noise, warmup, (0f0 , p_dict[:t_run]), p0)
#Run the solution to fruition
@time NetSol = solve(NetProb, SOSRI(), 
    abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
    save_idxs = [1:(Int64(p_dict[:nx]*p_dict[:ny]))...], 
    progress = true, progress_steps = 1
)

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