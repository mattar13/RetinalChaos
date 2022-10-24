#Load the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using RetinalChaos
using Plots
#using JLD2

#%% Run the simulations
#Step 1: Set up the network properties
print("[$(now())]: Setting up parameters, conditions, and network settings... ")
nx = ny = 125
net = RetinalChaos.T_PDE
#Step 2: Import the initial conditions
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))

#Step 3: Import the parameters
pars_dict = read_JSON("params/params.json")
p = pars_dict |> extract_dict

#Step 4: Determine the timespan
tspan = (0.0, 120e3)

#Step 5: Set up the problem
println("Complete")

#Step 6: Running the model
print("[$(now())]: Warming up the model for 60s... ")
prob = SDEProblem(net, noise, u0, (0.0, 60e3), p)
@time warmup = solve(prob, SOSRI(),
    abstol=2e-2, reltol=2e-2, maxiters=1e7,
    save_everystep=false, progress=true, progress_steps=1
);
println("Completed")

print("[$(Dates.now())]: Simulating up the model for $(round(tspan[end]/1000))s... ")
prob = SDEProblem(net, noise, warmup[end], tspan, p)
#Step 7: Run the model
@time NetSol = solve(prob, SROCK1(), dt=1.0,
    abstol=2e-2, reltol=0.2, maxiters=1e7,
    progress=true, progress_steps=1,
    save_idxs=[1:(nx*ny)...],
);
println("Completed")

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
gif(anim, "examples/animation.gif", fps=1000.0 / animate_dt)