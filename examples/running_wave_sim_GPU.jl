using Revise
using RetinalChaos
using Dates, Plots, JLD2

#Configure the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#Load the telegram client and env
dotenv("D:\\TelegramAccessEnv\\.env")
#Activate the GPU
RetinalChaos.CUDA.allowscalar(false)

#%% Load the needed files to run the model
save_file = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\$(Date(Dates.now()))_test1\\"
p_dict = read_JSON("params\\params.json", is_type = Dict{Symbol, Float32}) 
u_dict = read_JSON("params\\conds.json", is_type = Dict{Symbol, Float32})
NetSol = load_model(save_file, p_dict, u_dict, reset_model = true)
#%%

#%%Conduct the threshold analysis
thresholds = calculate_threshold(NetSol; dt = 500.0) #This takes really long
@save "$(save_file)\\thresholds.jld2" thresholds
BotNotify("{Wave} Finished running threshold analysis")

#%% Work in progress
timestamps = get_timestamps(NetSol, thresholds)


#%%
run_loop = false #so far set this section to ignore
if run_loop
    #%% Run multiple samples
    n = 4
    for val in LinRange(1.0, 20.0, 10)
        println(val)
        println("[$(Dates.now())]: Running sample ")
        # Setup the network simulation (Need to do this if accessing the saved file)
        print("[$(Dates.now())]: Setting up the model...")
        
        #Save the model, params, conditions, animations here
        save_file = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\$(Date(Dates.now()))_gCa_$val\\"
        if isdir(save_file) == false
            #The directory does not exist, we have to make it 
            mkdir(save_file)
        end
        p = read_JSON(params_file) 
        p[:g_Ca] = val
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
                save_idxs = [1:((p[:nx]*p[:ny])|>Int64)...], 
                progress = true, progress_steps = 1
            )
        println("Completed")

        # Save the solution, must be on drive first
        print("[$(Dates.now())]: Saving the simulation...")
        JLD2.@save "$(save_file)\\sol.jld2" NetSol
        
        # Plotting animation
        anim = @animate for t = 1.0:60.0:NetSol.t[end]
            println("Animating frame $t")
            frame_i = reshape(NetSol(t) |> Array, (p[:nx]|>Int64, p[:ny]|>Int64))
            heatmap(frame_i, ratio = :equal, grid = false,
                    xaxis = "", yaxis = "", xlims = (0, p[:nx]|>Int64), ylims = (0, p[:ny]|>Int64),
                    c = :curl, clims = (-70.0, 0.0),
            )
        end
        #
        gif(anim, "$(save_file)\\animation.gif", fps = 40.0)
    end
end