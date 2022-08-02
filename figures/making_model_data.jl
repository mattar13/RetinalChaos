using Revise
using RetinalChaos
using Plots
include("figure_setup.jl")
# Run 3 models
#%% Model 1: Regular Baseline model 
print("[$(now())]: Setting up parameters, conditions, and network settings... ")
nx = ny = 64
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
pars_dict = read_JSON("params/params.json")
p = pars_dict |> extract_dict
tspan = (0.0, 120e3)
println("Complete")
print("[$(now())]: Warming up the model for 60s... ")
prob = SDEProblem(T_PDE, noise, u0, (0.0, 60e3), p)
@time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1);
println("Completed")
print("[$(now())]: Simulating up the model for $(round(tspan[end]/1000))s... ")
prob = SDEProblem(T_PDE, noise, warmup[end], tspan, p)
@time sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...]);
println("Completed")
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\wave_model"
if !isdir(loc) #If the directory doesn't exist, make it
     println("directory doesn't exist. Making it")
     mkdir(loc)
end
animate_dt = 60.0
anim = @animate for t = 1.0:animate_dt:sol.t[end]
     println("[$(now())]: Animating simulation $(t) out of $(sol.t[end])...")
     frame_i = reshape(sol(t) |> Array, (nx, ny))
     heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
end
gif(anim, "$(loc)/regular_animation.gif", fps=1000.0 / animate_dt)
timestamps, data = timeseries_analysis(sol, loc)
hist_plot = plot_histograms(data, loc)

#%% Model 2: Blocked Neurotransmission 
print("[$(now())]: Setting up parameters, conditions, and network settings... ")
nx = ny = 64
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
pars_dict = read_JSON("params/params.json")
pars_dict[:g_ACh] = 0.0 # Block all Acetylcholine receptors
pars_dict[:g_GABA] = 0.0 #Block all GABA receptors
p = pars_dict |> extract_dict
tspan = (0.0, 120e3)
println("Complete")
print("[$(now())]: Warming up the model for 60s... ")
prob = SDEProblem(T_PDE, noise, u0, (0.0, 60e3), p)
@time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1);
println("Completed")
print("[$(now())]: Simulating up the model for $(round(tspan[end]/1000))s... ")
prob = SDEProblem(T_PDE, noise, warmup[end], tspan, p)
@time sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...]);
println("Completed")
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\isolated_model"
if !isdir(loc) #If the directory doesn't exist, make it
     println("directory doesn't exist. Making it")
     mkdir(loc)
end
animate_dt = 60.0
anim = @animate for t = 1.0:animate_dt:sol.t[end]
     println("[$(now())]: Animating simulation $(t) out of $(sol.t[end])...")
     frame_i = reshape(sol(t) |> Array, (nx, ny))
     heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
end
gif(anim, "$(loc)/regular_animation.gif", fps=1000.0 / animate_dt)
timestamps, data = timeseries_analysis(sol, loc)
hist_plot = plot_histograms(data, loc)

#%% Model 3: No GABA
print("[$(now())]: Setting up parameters, conditions, and network settings... ")
nx = ny = 64
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
pars_dict = read_JSON("params/params.json")
pars_dict[:g_GABA] = 0.0 #Block all GABA receptors
p = pars_dict |> extract_dict
tspan = (0.0, 120e3)
println("Complete")
print("[$(now())]: Warming up the model for 60s... ")
prob = SDEProblem(T_PDE, noise, u0, (0.0, 120e3), p) #Extend the warmup phase to get more burst behavior
@time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1);
println("Completed")
print("[$(now())]: Simulating up the model for $(round(tspan[end]/1000))s... ")
prob = SDEProblem(T_PDE, noise, warmup[end], tspan, p)
@time sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...]);
println("Completed")
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\no_GABA_model"
if !isdir(loc) #If the directory doesn't exist, make it
     println("directory doesn't exist. Making it")
     mkdir(loc)
end
animate_dt = 60.0
anim = @animate for t = 1.0:animate_dt:sol.t[end]
     println("[$(now())]: Animating simulation $(t) out of $(sol.t[end])...")
     frame_i = reshape(sol(t) |> Array, (nx, ny))
     heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
end

gif(anim, "$(loc)/regular_animation.gif", fps=1000.0 / animate_dt)
timestamps, data = timeseries_analysis(sol, loc)
hist_plot = plot_histograms(data, loc)
data["SpikeDurs"]
#%% lets make a loop that will change gCa but also lack GABA
#gCa Spiking range = [7.0-9.0]
param = :g_K
for val in LinRange(4.0, 8.0, 4) #This is the range of values
     print("[$(now())]: Setting up parameters, conditions, and network settings... ")
     nx = ny = 64
     conds_dict = read_JSON("params/conds.json")
     u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
     pars_dict = read_JSON("params/params.json")
     pars_dict[param] = val
     pars_dict[:g_GABA] = 0.0 #Inhibit GABA
     p = pars_dict |> extract_dict
     tspan = (0.0, 120e3)
     println("Complete")
     print("[$(now())]: Warming up the model for 60s... ")
     prob = SDEProblem(T_PDE, noise, u0, (0.0, 60e3), p)
     @time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1)
     println("Completed")
     print("[$(now())]: Simulating up the model for $(round(tspan[end]/1000))s... ")
     prob = SDEProblem(T_PDE, noise, warmup[end], tspan, p)
     @time sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...])
     println("Completed")
     name = "GABA_0_$(param)_$(round(val, digits = 2))"
     loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/experiments/$(name)"
     if !isdir(loc) #If the directory doesn't exist, make it
          println("directory doesn't exist. Making it")
          mkdir(loc)
     end
     animate_dt = 60.0
     anim = @animate for t = 1.0:animate_dt:sol.t[end]
          println("[$(now())]: Animating simulation $(t) out of $(sol.t[end])...")
          frame_i = reshape(sol(t) |> Array, (nx, ny))
          heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
     end
     title = "$(param)_$(round(val, digits=2))"
     gif(anim, "$(loc)/$(title)_animation.gif", fps=1000.0 / animate_dt)
     timestamps, data = timeseries_analysis(sol, loc; tstamps_name="$(title)_timestamps", data_name="$(title)_data")
     hist_plot = plot_histograms(data, loc; name="$(title)_histogram")
end