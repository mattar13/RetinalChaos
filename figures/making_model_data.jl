using Revise
using RetinalChaos
using Plots
import RetinalChaos: write_JSON, plot_histograms
import Plots: @animate

# Run 2 models
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
@save "$(loc)/regular_sol.jld2" sol
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
@save "$(loc)/null_sol.jld2" sol
animate_dt = 60.0
anim = @animate for t = 1.0:animate_dt:sol.t[end]
     println("[$(now())]: Animating simulation $(t) out of $(sol.t[end])...")
     frame_i = reshape(sol(t) |> Array, (nx, ny))
     heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
end
gif(anim, "$(loc)/regular_animation.gif", fps=1000.0 / animate_dt)
timestamps, data = timeseries_analysis(sol, loc)
hist_plot = plot_histograms(data, loc)

#%% make a for loop that goes through some parameters (we can plot)
param = :g_Ca
for val in LinRange(5.0, 10.0, 10) #This is the range of values
     print("[$(now())]: Setting up parameters, conditions, and network settings... ")
     nx = ny = 50
     conds_dict = read_JSON("params/conds.json")
     u0 = extract_dict(conds_dict, t_conds, dims=(nx, ny))
     pars_dict = read_JSON("params/params.json")
     pars_dict[param] = val
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
     loc = raw"C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data/range_gCa/"
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