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
net = GABA_PDE
#Step 2: Import the initial conditions
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))

#Step 3: Import the parameters
pars_dict = read_JSON("params/params.json")
p = pars_dict |> extract_dict

#Step 4: Determine the timespan
tspan = (0.0, 120e3)

#Step 5: Set up the problem
prob = SDEProblem(net, noise, u0, tspan, p)

#Step 6: Running the model
@time warmup = solve(prob, SOSRI(),
    abstol=2e-2, reltol=2e-2, maxiters=1e7,
    save_everystep=false, progress=true, progress_steps=1
)
tspan = (0.0, 120e3)
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
        c=:curl, clims=(-90.0, 0.0),
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