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
pars_dict[:C_m] = 13.6
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
pars_dict[:C_m] = 13.6
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

#%% Model 4 ECl Differential
print("[$(now())]: Setting up parameters, conditions, and network settings... ")
nx = ny = 64
conds_dict = read_JSON("params/conds.json")
u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
pars_dict = read_JSON("params/params.json")
pars_dict[:E_GABA] = -55.0 #Block all GABA receptors
p = pars_dict |> extract_dict
tspan = (0.0, 120e3)
println("Complete")
print("[$(now())]: Warming up the model for 60s... ")
prob = SDEProblem(T_PDE, RetinalChaos.noise, u0, (0.0, 120e3), p) #Extend the warmup phase to get more burst behavior
@time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1);
println("Completed")
print("[$(now())]: Simulating up the model for $(round(tspan[end]/1000))s... ")
prob = SDEProblem(T_PDE, RetinalChaos.noise, warmup[end], tspan, p)
@time sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...]);
println("Completed")
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\ECl55_model"
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

#%% if we just wanted to open the data to plot a heatmap figure1_ModelVariables
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
nx = ny = 64
animate_dt = 60.0

data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
isolated_path = "$(data_root)/isolated_model"
noGABA_path = "$(data_root)/no_GABA_model"
wave_path = "$(data_root)/wave_model"
ECl55_path = "$(data_root)/ECl55_model"

dataISO = load("$(isolated_path)/data.jld2")
dataNG = load("$(noGABA_path)/data.jld2")
dataWAVE = load("$(wave_path)/data.jld2")
dataECl = load("$(ECl55_path)/data.jld2")

solISO = reshape(dataISO["DataArray"], (nx, ny, size(dataISO["DataArray"], 2)))
solNG = reshape(dataNG["DataArray"], (nx, ny, size(dataNG["DataArray"], 2)))
solWAVE = reshape(dataWAVE["DataArray"], (nx, ny, size(dataWAVE["DataArray"], 2)))
solECl = reshape(dataECl["DataArray"], (nx, ny, size(dataECl["DataArray"], 2)))

save_loc = "C:\\Users\\mtarc\\The University of Akron\\Renna Lab - General\\Journal Submissions\\2022 A Computational Model - Sci. Rep\\Submission 1\\Figures\\"
#%% Save the isolated data model
anim = @animate for t = 1.0:animate_dt:size(solISO, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solISO, 2))...")
     frame_i = solISO[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S2 Neurotransmission Blocked Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the no GABA simulation
anim = @animate for t = 1.0:animate_dt:size(solNG, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solNG, 2))...")
     frame_i = solNG[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S3 NoGABA Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the wave simulation
anim = @animate for t = 1.0:animate_dt:size(solWAVE, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solNG, 2))...")
     frame_i = solWAVE[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S4 Wave Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the ECl-55 simulation
anim = @animate for t = 1.0:animate_dt:size(solECl, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solECl, 2))...")
     frame_i = solECl[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S5 ECl_-55 Simulation.gif", fps=1000.0 / animate_dt)