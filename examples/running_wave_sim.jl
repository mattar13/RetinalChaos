using RetinalChaos
using Dates
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")
save_path = "C:\\Users\\RennaLabSA1\\Documents\\modelling"

open(joinpath(save_path, "log.txt"), "w") do logf
    
println(logf, "[$(Dates.now())]: File set up")
#%%
println(logf, "[$(Dates.now())]: Beginning simulion")
println("[$(Dates.now())]: Beginning simulion")
for sig in range(0.01, 0.25, length = 10), ga in range(0.3, 2.0, length = 10)

    println(logf, "[$(Dates.now())]:Running test for gACh = $(ga) sigma = $(sig)")
    println("[$(Dates.now())]:Running test for gACh = $(ga) sigma = $(sig)")
    nx = ny = 125; 
    p = read_JSON(params_file);
    p[:σ] = sig
    p[:g_ACh] = ga
    p = extract_dict(p, tar_pars);
    u0_network = extract_dict(read_JSON(conds_file), tar_conds, (nx, ny));
    net = Network(nx, ny; μ = 0.60, version = :gACh) #μ is the probability a cell is capable of being active
    var = 1 #This is the variable we are interested in
    save_idxs = [var*1:var*(nx*ny)...] #This is a list of indexes of all the variables we want to plot

    #%% Lets warm up the solution first
    PDEprob = SDEProblem(net, noise, u0_network, (0.0, 10e3), p)
    print(logf, "[$(Dates.now())]: Warming up solution:")
    print("[$(Dates.now())]: Warming up solution:")
    @time PDEsol = solve(PDEprob, SOSRI(), 
        abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
        save_everystep = false, progress = true, progress_steps = 1
        )
    # extract the last solution and use it as the initial condition for future models 
    # This may take a long time to run
    print(logf, "[$(Dates.now())]: Running model:")
    print("[$(Dates.now())]: Running model:")
    ui_network = PDEsol[end]
    PDEprob = SDEProblem(net, noise, ui_network, (0.0, 300e3), p)
    @time PDEsol = solve(PDEprob, SOSRI(), 
        abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
        save_idxs = save_idxs, progress = true, progress_steps = 1
        )
    #%% Plot as a .gif
    save_wave = joinpath(save_path, "sig_$(sig)_gACh_$(ga)_wave_plot.gif")
    save_plot = joinpath(save_path, "sig_$(sig)_gACh_$(ga)_2D_plot.gif")
    
    println(logf, "[$(Dates.now())]: Plotting results")
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

    #%%
    plt = plot()
    for i in 1:20
        println(i)
        plot!(plt, PDEsol.t, PDEsol[i,:])
    end 
    savefig(plt, save_plot)
end
end