using RetinalChaos
using Dates, Plots, CUDA
using JLD2
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")
#Load the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
CUDA.allowscalar(false)

#%% Setup the network simulation (Need to do this if accessing the saved file)
print("[$(Dates.now())]: Setting up the model...")
#Save the model, params, conditions, animations here
save_file = "data\\$(Date(Dates.now()))_gACh_mu_0.15\\"
if isdir(save_file) == false
    #The directory does not exist, we have to make it 
    mkdir(save_file)
end
p = read_JSON(params_file) 
u0 = read_JSON(conds_file);
#Save the params and ics to load later. 
write_JSON(p, "$(save_file)\\params.json") #write the parameters to save for later
write_JSON(u0, "$(save_file)\\conds.json") #write the ics to save for later
net = Network(p[:nx]|>Int64, p[:ny]|>Int64; μ = p[:μ], version = :gACh, gpu = true) #This branch uses GPU mode
p_net = extract_dict(p);
u0_net = extract_dict(u0, p[:nx]|>Int64, p[:ny]|>Int64) |> cu;
warmup = (0.0|>Float32, p[:t_warm]|>Float32)
tspan  = (0.0|>Float32, p[:t_run]|>Float32)
NetProb = SDEProblem(net, noise, u0_net, warmup, p_net)
println("Completed")


#%% Lets warm up the solution first (using GPU if available)
print("[$(Dates.now())]: Warming up the solution: ")
@time NetSol = solve(NetProb, SOSRI(), 
        abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
        progress = true, progress_steps = 1, 
        save_everystep = false
    )
warmup_ics = NetSol[end]
println("Completed")

#%% Save or load the warmed up solution
print("[$(Dates.now())]: Loading or saving solution...")
JLD2.@save "$(save_file)\\warmup_ics.jld2" warmup_ics
println("Completed")

#%% Run the simulation
print("[$(Dates.now())]: Running the simulation... ")
NetProb = SDEProblem(net, noise, warmup_ics, tspan, p_net)
@time NetSol = solve(NetProb, SOSRI(), 
        abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
        save_idxs = [1:(nx*ny)...], 
        progress = true, progress_steps = 1
    )
println("Completed")

#%% Save the solution, must be on drive first
print("[$(Dates.now())]: Saving the simulation...")
JLD2.@save "$(save_file)\\sol.jld2" NetSol
#%% Plotting animation
anim = @animate for t = 1.0:60.0:NetSol.t[end]
    println("Animating frame $t")
    frame_i = reshape(NetSol(t) |> Array, (nx, ny))
    heatmap(frame_i, ratio = :equal, grid = false,
            xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
            c = :curl, clims = (-70.0, 0.0),
    )
end
#%%
gif(anim, "$(save_file)\\animation.gif", fps = 40.0)

#%%We can load the data here 
file_root = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\2021-04-28_gACh_mu_0.2"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

print("[$(Dates.now())]: Loading the model...")
p = read_JSON(params_file) 
#Set up the initial conditions
u0 = read_JSON(conds_file);
net = Network(p[:nx], p[:ny]; μ = p[:μ], version = :gACh, gpu = true) #This branch uses GPU mode
p_net = extract_dict(p);
u0_net = extract_dict(u0, p[:nx], p[:ny]) |> cu;
warmup = (0.0|>Float32, 300e3|>Float32)
tspan  = (0.0|>Float32 , 60e3|>Float32)
NetProb = SDEProblem(net, noise, u0_net, warmup, p_net)
JLD2.@load "$(file_root)\\sol.jld2" NetSol
println("Completed")

#%%
save_arr = Array(NetSol)
#%%
thresholds = RetinalChaos.calculate_threshold(NetSol) #This takes really long
#%%
@save "$(save_file)\\thresholds.jld2" thresholds
#%%
ts = RetinalChaos.timescale_analysis(NetSol, thresholds)
@save "$(save_file)\\ts_analysis.jld2"


#%% Run multiple samples
n = 4
for sample in 1:n
    println("[$(Dates.now())]: Running sample ")
    # Setup the network simulation (Need to do this if accessing the saved file)
    print("[$(Dates.now())]: Setting up the model...")
    
    #Save the model, params, conditions, animations here
    save_file = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\$(Date(Dates.now()))_sample_$sample\\"
    if isdir(save_file) == false
        #The directory does not exist, we have to make it 
        mkdir(save_file)
    end
    p = read_JSON(params_file) 
    u0 = read_JSON(conds_file);

    #Save the params and ics to load later. 
    write_JSON(p, "$(save_file)\\params.json") #write the parameters to save for later
    write_JSON(u0, "$(save_file)\\conds.json") #write the ics to save for later
    net = Network(p[:nx]|>Int64, p[:ny]|>Int64; μ = p[:μ], version = :gACh, gpu = true) #This branch uses GPU mode
    p_net = extract_dict(p);
    u0_net = extract_dict(u0, p[:nx]|>Int64, p[:ny]|>Int64) |> cu;
    warmup = (0.0|>Float32, p[:t_warm]|>Float32)
    tspan  = (0.0|>Float32, p[:t_run]|>Float32)
    NetProb = SDEProblem(net, noise, u0_net, warmup, p_net)
    println("Completed")

    # Lets warm up the solution first (using GPU if available)
    print("[$(Dates.now())]: Warming up the solution: ")
    @time NetSol = solve(NetProb, SOSRI(), 
            abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
            progress = true, progress_steps = 1, 
            save_everystep = false
        )
    warmup_ics = NetSol[end]
    println("Completed")

    # Save or load the warmed up solution
    print("[$(Dates.now())]: Loading or saving solution...")
    JLD2.@save "$(save_file)\\warmup_ics.jld2" warmup_ics
    println("Completed")

    # Run the simulation
    print("[$(Dates.now())]: Running the simulation... ")
    NetProb = SDEProblem(net, noise, warmup_ics, tspan, p_net)
    @time NetSol = solve(NetProb, SOSRI(), 
            abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
            save_idxs = [1:(nx*ny)...], 
            progress = true, progress_steps = 1
        )
    println("Completed")

    # Save the solution, must be on drive first
    print("[$(Dates.now())]: Saving the simulation...")
    JLD2.@save "$(save_file)\\sol.jld2" NetSol
    
    # Plotting animation
    anim = @animate for t = 1.0:60.0:NetSol.t[end]
        println("Animating frame $t")
        frame_i = reshape(NetSol(t) |> Array, (nx, ny))
        heatmap(frame_i, ratio = :equal, grid = false,
                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
                c = :curl, clims = (-70.0, 0.0),
        )
    end
    #
    gif(anim, "$(save_file)\\animation.gif", fps = 40.0)
end