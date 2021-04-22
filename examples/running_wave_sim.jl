using RetinalChaos
using Dates, Plots
using JLD2
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")


#%% Run a network simulation and save it
nx = 64 
ny = 64; 
p = read_JSON(params_file) 
p[:σ] = 0.1
p[:τw] = 800.0
#Set up the initial conditions
u0 = read_JSON(conds_file);
net = Network(nx, ny; μ = 0.15, version = :ρ) 
p_net = extract_dict(p, tar_pars);
u0_net = extract_dict(u0, tar_conds, (nx, ny)) |> cu;
NetProb = SDEProblem(net, noise, u0_net, (0.0, 60e3), p_net)
#%%
#Lets warm up the solution first
@time NetSol = solve(NetProb, SOSRI(), 
    abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
    save_everystep = false)

#Run the simulation
NetProb = SDEProblem(net, noise, NetSol[end], (0.0, 120e3), p_net)
@time NetSol = solve(NetProb, SOSRI(), 
    abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
    save_idxs = [1:(nx*ny)...], 
    progress = true, progress_steps = 1
    )

#%% Save the solution, must be on drive first
JLD2.@save "$(Date(Dates.now()))_sol.jld2" NetSol
#%% Plotting stuff
anim = @animate for t = 1.0:50.0:PDEsol.t[end]
    println("Animating frame $t")
    frame_i = reshape(PDEsol(t), (nx, ny))
    heatmap(frame_i, ratio = :equal, grid = false,
            xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
            c = :curl, clims = (-70.0, 0.0),
    )
end
gif(anim, save_wave, fps = 20)

#%% Extracting wave data from the simulation
open_dir = "C:\\Users\\mtarc\\OneDrive\\Documents\\GithubRepositories\\RetinalChaos\\figures\\2021-04-19_sol.jld2"
JLD2.@load open_dir NetSol
#%%
thresholds = RetinalChaos.calculate_threshold(NetSol) #This takes really long
@save "figures\\thresholds.jld2" thresholds
#%%
ts = RetinalChaos.timescale_analysis(NetSol, thresholds)
@save "figures\\ts_analysis.jld2"

#%% Running a multiple sim on the Lab computer
save_path = "C:\\Users\\RennaLabSA1\\Documents\\modelling"
println("[$(Dates.now())]: Beginning simulion")
ga = 1.1
sig = 0.25
for rho in range(0.1, 1.0, length = 50)#for sig in range(0.01, 1.0, length = 10)#, ga in range(0.1, 1.1, length = 2)
    println("[$(Dates.now())]:Running test for gACh = $(ga) sigma = $(sig)")
    nx = ny = 64; 
    p = read_JSON(params_file);
    p[:σ] = sig
    p[:g_ACh] = ga
    p = extract_dict(p, tar_pars);
    u0_network = extract_dict(read_JSON(conds_file), tar_conds, (nx, ny));
    net = Network(nx, ny; μ = rho, version = :ρ) #μ is the probability a cell is capable of being active
    var = 1 #This is the variable we are interested in
    save_idxs = [var*1:var*(nx*ny)...] #This is a list of indexes of all the variables we want to plot

    # Lets warm up the solution first
    PDEprob = SDEProblem(net, noise, u0_network, (0.0, 60e3), p)
    print("[$(Dates.now())]: Warming up solution:")
    @time PDEsol = solve(PDEprob, SOSRI(), 
        abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
        save_everystep = false, progress = true, progress_steps = 1
        )
    # extract the last solution and use it as the initial condition for future models 
    # This may take a long time to run
    print("[$(Dates.now())]: Running model:")
    ui_network = PDEsol[end]
    PDEprob = SDEProblem(net, noise, ui_network, (0.0, 120e3), p)
    @time PDEsol = solve(PDEprob, SOSRI(), 
        abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
        save_idxs = save_idxs, progress = true, progress_steps = 1
        )
    # Plot as a .gif
    save_wave = joinpath(save_path, "sig_$(replace(string(sig), "." => "_"))_gACh_$(replace(string(ga), "." => "_"))_rho_$(replace(string(rho), "." => "_"))wave_plot.gif")
    save_plot = joinpath(save_path, "sig_$(replace(string(sig), "." => "_"))_gACh_$(replace(string(ga), "." => "_"))_rho_$(replace(string(rho), "." => "_"))2D_plot.png")
    
    println("[$(Dates.now())]: Plotting results")
    dFrame = 50.0
    t_rng = collect(1.0:dFrame:PDEsol.t[end])
    anim = @animate for t = t_rng
        println("Animating frame $t")
        frame_i = reshape(PDEsol(t), (nx, ny))
        heatmap(frame_i, ratio = :equal, grid = false,
                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
                c = :curl, clims = (-70.0, 0.0),
        )
    end
    gif(anim, save_wave, fps = 20)

    #
    plt = plot()
    for i in 1:20
        println(i)
        plot!(plt, PDEsol.t, PDEsol[i,:])
    end 
    savefig(plt, save_plot)
end