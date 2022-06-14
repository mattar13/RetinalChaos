#Load the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using RetinalChaos
using Plots
#using JLD2

#%% Run the simulations
#Step 1: Set up the network properties
nx = ny = 50

#Step 1b: If using binomial nullification set that up now
b = Binomial(1, 0.90)
null = Array{Float64}(rand(b, nx, ny))
heatmap(null; ratio=:equal)
net = (dU, U, p, t) -> GABA_PDE(dU, U, p, t)

#Step 2: Import the initial conditions
conds_dict = read_JSON("params/GABA_conds.json")
u0 = extract_dict(conds_dict, GABA_conds, dims=(nx, ny))

#Step 3: Import the parameters
pars_dict = read_JSON("params/GABA_params.json")
pars_dict[:g_ACh] = 0.215
pars_dict[:g_GABA] = 1.1
p = extract_dict(pars_dict, GABA_pars)

#Step 4: Determine the timespan
tspan = (0.0, 120e3)
#Step 5: Set up the problem
prob = SDEProblem(net, noise, u0, tspan, p)


#Step 6: Running the model
@time warmup = solve(prob, SOSRI(),
    abstol=2e-2, reltol=2e-2, maxiters=1e7,
    save_everystep=false, progress=true, progress_steps=1
)
prob = SDEProblem(net, noise, warmup[end], tspan, p)

#Step 7: Run the model
@time NetSol = solve(prob, SROCK1(), dt=1.0,
    abstol=2e-2, reltol=0.2, maxiters=1e7,
    progress=true, progress_steps=1,
    save_idxs=[1:(nx*ny)...],
)

# Step 7: Animate the solution
animate_dt = 60.0
anim = @animate for t = 1.0:animate_dt:NetSol.t[end]
    println("[$(now())]: Animating simulation...")
    frame_i = reshape(NetSol(t) |> Array, (nx, ny))
    heatmap(frame_i, ratio=:equal, grid=false,
        xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
        c=:curl, clims=(-70.0, 0.0),
    )
end
gif(anim, "animation.gif", fps=1000.0 / animate_dt)
#%% Can we plot a solution
size(NetSol)
plot(NetSol.t, t -> NetSol(t, idxs=[1, 2, 3])')
#%% Save or load the warmed up solution
print("[$(Dates.now())]: Loading or saving solution...")
JLD2.@save "$(save_file)\\warmup_ics.jld2" warmup_ics
#JLD2.@load "$(save_file)\\warmup_ics.jld2" warmup_ics
println("Completed")

#%% Run the simulation
print("[$(Dates.now())]: Running the simulation... ")
NetProb = SDEProblem(net, noise, warmup_ics, tspan, p_net)
@time NetSol = solve(NetProb, SOSRI(),
    abstol=2e-2, reltol=0.2, maxiters=1e7,
    save_idxs=[1:(nx*ny)...],
    progress=true, progress_steps=1
)
println("Completed")

#%% Save the solution, must be on drive first
print("[$(Dates.now())]: Saving the simulation...")
JLD2.@save "$(save_file)\\sol.jld2" NetSol
println("Completed")
#%%
plot(NetSol, idxs=1)
#%% Plotting animation
anim = @animate for t = 1.0:10.0:NetSol.t[end]
    println("Animating frame $t")
    frame_i = reshape(NetSol(t) |> Array, (nx, ny))
    heatmap(frame_i, ratio=:equal, grid=false,
        xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
        c=:curl, clims=(-70.0, 0.0),
    )
end
gif(anim, "$(save_file)\\animation.gif", fps=20)

#%%
thresholds = RetinalChaos.calculate_threshold(NetSol) #This takes really long
@save "$(save_file)\\thresholds.jld2" thresholds
#%%
ts = RetinalChaos.timescale_analysis(NetSol, thresholds)
@save "$(save_file)\\ts_analysis.jld2"

#= %% Running a multiple sim on the Lab computer
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
end =#