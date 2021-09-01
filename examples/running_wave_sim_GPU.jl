using Revise
using RetinalChaos
using Dates, Plots, BSON
using StatsBase, StatsPlots

#Configure the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
import RetinalChaos.param_path

#Load the telegram client and env
dotenv("D:\\TelegramAccessEnv\\.env")
#Activate the GPU
RetinalChaos.CUDA.allowscalar(false)

#%% Run a simulation on it's own
param_test = :μ #Set a single parameter
x = 0.50 #Set it's value
NetSol = nothing #this allows the solution to be used outside of the try loop
save_path = nothing
try
    BotNotify("{Waves} Running simulation for $(param_test) = $x started")
    u_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/conds.json") #Load initial conditions
    p_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/params.json") #Load the parameters
    ################# If you want to modify parameters do it here ###############
    p_dict[:t_run] = 300e3 #We want to make the run time longer
    p_dict[param_test] = x
    ################# If you want to modify parameters do it here ###############
    save_path = "E:\\Data\\Modelling\\$(param_test)_$(round(Int64, x*100))" #Specify the save path
    abs_tol = 2e-2
    rel_tol = 2e-2 #try changing this to get less instability
    NetSol = run_model(save_path, p_dict, u_dict; abstol = abs_tol, reltol = rel_tol) #This is a succinct function for running the model
    BotNotify("{Waves} Running simulation for $(param_test) = $x completed")
    #This is how we can calculate how much space the data will take up
catch error
    BotNotify("{Waves} An error has occurred $(typeof(error))")
    throw(error)
end

#%% Save the solution
try
    save_solution(NetSol, save_path, partitions = 50) #Save the solution in 50 seperate parts
    BotNotify("{Waves} Saving simulation for $(param_test) = $x completed")
catch error
    BotNotify("{Waves} Saving simulation for $(param_test) = $x failed")
    throw(error)
end
#%%
BotNotify("{Waves} Saving simulation for $(param_test) = $x completed")
sol = load_solution(save_path) #Reload the solution for analysis

#BotNotify("{Waves} Animating simulation for $(param_test) = $x completed")
timeseries_analysis(sol, save_path, dt = 1.0, verbose = true) #conduct a timeseries analysis
animate_solution(sol, save_path) #Animate the solution
BotNotify("{Waves} Timeseries analysis completed")
data = BSON.load("$(save_path)\\data.bson")

#%% Run a longer simulation loop
#param_range = LinRange(0.05, 1.0, 25) #We want to rerun this exp with wave extraction
param_test = :μ
param_range = [0.0, 0.125, 0.25, 0.50]
for x in param_range
    println("Running simulation for $(param_test) = $x")
    try
        BotNotify("{Waves} Running simulation for $(string(param_test)) = $x started")
        u_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/conds.json")
        p_dict = read_JSON(Dict{Symbol, Float32}, "$(param_path)/params.json")
        ############### If you want to modify parameters do it here ###############
        p_dict[:t_run] = 300e3 #We want to make the run time longer
        p_dict[param_test] = x
        ############### If you want to modify parameters do it here ###############
        save_path = "E:\\Data\\Modelling\\mu_experiment\\$(string(param_test))_$(round(Int64, x*100))\\"
        NetSol = run_model(save_path, p_dict, u_dict)
        BotNotify("{Waves} Running simulation for $(param_test) = $x completed")
        save_solution(NetSol, save_path, partitions = 10)
        BotNotify("{Waves} Saving simulation for $(param_test) = $x completed")
        sol = load_solution(save_path)
        animate_solution(sol, save_path)
        BotNotify("{Waves} Animating simulation for $(param_test) = $x completed")
        timeseries_analysis(sol, save_path)
        BotNotify("{Waves} Timeseries analysis for $(param_test) = $x completed")
    catch error
        println(error)
        BotNotify("{Waves} An error has occurred $(typeof(error))")
    end
end