using Revise
using RetinalChaos
using DataFrames
include("figure_setup.jl")

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using ProgressMeter
# Run 3 models

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling"
#%% Model 1: Regular Baseline model 
loc = "$(data_root)\\wave_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
@time sol = run_model(p_dict, u_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7);

#%% Model 2: Blocked Neurotransmission 
loc = "$(data_root)\\isolated_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:g_ACh] = 0.0 # Block all Acetylcholine receptors
p_dict[:g_GABA] = 0.0 #Block all GABA receptors
@time sol = run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Model 3: No GABA
loc = "$(data_root)\\no_GABA_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:g_GABA] = 0.0 #Block all GABA receptors
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Model 4 ECl Differential
loc = "$(data_root)\\ECl55_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:E_Cl] = -55.0 #Block all GABA receptors
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Model 5 ECl more hyperpolarizing Differential
loc = "$(data_root)\\ECl75_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:E_Cl] = -75.0 #Block all GABA receptors
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Plot the wave simulation
nx = ny = 64
nt = 120001
wave_path = "$(data_root)/FastI_EDiffusion"
dataWAVE = load("$(wave_path)/data.jld2")
timestampsWAVE = load("$(wave_path)/timestamps.jld2")
burstWAVE = RetinalChaos.extract_burstmap(nx,ny,nt,timestampsWAVE["Bursts"])
WaveSegmentation(burstWAVE, wave_path)

isolated_path = "$(data_root)/isolated_model"
dataISO = load("$(isolated_path)/data.jld2")
timestampsISO = load("$(isolated_path)/timestamps.jld2")
burstISO = RetinalChaos.extract_burstmap(nx, ny, nt, timestampsISO["Bursts"])
WaveSegmentation(burstISO, isolated_path)

noGABA_path = "$(data_root)/no_GABA_model"
dataNG = load("$(noGABA_path)/data.jld2")
timestampsNG = load("$(noGABA_path)/timestamps.jld2")
burstNG = RetinalChaos.extract_burstmap(nx, ny, nt, timestampsNG["Bursts"])
WaveSegmentation(burstNG, noGABA_path)

ECl55_path = "$(data_root)/ECl55_model"
dataECl = load("$(ECl55_path)/data.jld2")
timestampsECl = load("$(noGABA_path)/timestamps.jld2")
burstECl = RetinalChaos.extract_burstmap(nx, ny, nt, timestampsECl["Bursts"])
WaveSegmentation(burstECl, ECl55_path)

#%% Animations Wave
solWAVE = reshape(dataWAVE["DataArray"], (nx, ny, size(dataWAVE["DataArray"], 2)))
anim = @animate for t = 1.0:animate_dt:size(solWAVE, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solWAVE, 3))...")
     frame_i = solWAVE[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S6 Wave Simulation.gif", fps=1000.0 / animate_dt)

#%% Animations isolated path
solISO = reshape(dataISO["DataArray"], (nx, ny, size(dataISO["DataArray"], 2)))
anim = @animate for t = 1.0:animate_dt:size(solISO, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solISO, 2))...")
     frame_i = solISO[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S4 Neurotransmission Blocked Simulation.gif", fps=1000.0 / animate_dt)

#%% No GABA path
solNG = reshape(dataNG["DataArray"], (nx, ny, size(dataNG["DataArray"], 2)))
anim = @animate for t = 1.0:animate_dt:size(solNG, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solNG, 2))...")
     frame_i = solNG[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S5 NoGABA Simulation.gif", fps=1000.0 / animate_dt)


#%% Plot the ECl-55 simulation
solECl = reshape(dataECl["DataArray"], (nx, ny, size(dataECl["DataArray"], 2)))
anim = @animate for t = 1.0:animate_dt:size(solECl, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solECl, 2))...")
     frame_i = solECl[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S7 ECl -55 Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the ECl-75 simulation
solECl = reshape(dataECl["DataArray"], (nx, ny, size(dataECl["DataArray"], 2)))
anim = @animate for t = 1.0:animate_dt:size(solECl, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solECl, 2))...")
     frame_i = solECl[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S7 ECl -55 Simulation.gif", fps=1000.0 / animate_dt)


#Open up data for the Diffusion experiments==================================================================================================================#

#%% What are the effects of changing the synaptic release? 
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\FastEdiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:De] = 0.01 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#Fast inhibitory ==================================================================================================================#
#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\FastIDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:Di] = 0.01 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\FastI_EDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:De] = 0.01 #Make synaptic transmission faster
p_dict[:Di] = 0.01 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% What are the effects of changing the synaptic release? 
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\SlowEdiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:De] = 0.001 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\SlowIDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:Di] = 0.001 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\SlowI_EDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:Di] = 0.001 #Make synaptic transmission faster
p_dict[:De] = 0.001 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc, alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Save the isolated data model
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
FASTe_path = "$(data_root)/isolated_model"
FASTi_path = "$(data_root)/no_GABA_model"
FASTei_path = "$(data_root)/wave_model"


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

anim = @animate for t = 1.0:animate_dt:size(solISO, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solISO, 2))...")
     frame_i = solISO[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S3 Neurotransmission Blocked Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the no GABA simulation
anim = @animate for t = 1.0:animate_dt:size(solNG, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solNG, 2))...")
     frame_i = solNG[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S4 NoGABA Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the wave simulation
anim = @animate for t = 1.0:animate_dt:size(solWAVE, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solNG, 2))...")
     frame_i = solWAVE[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S5 Wave Simulation.gif", fps=1000.0 / animate_dt)

#%% Plot the ECl-55 simulation
anim = @animate for t = 1.0:animate_dt:size(solECl, 3)
     println("[$(now())]: Animating simulation $(t) out of $(size(solECl, 2))...")
     frame_i = solECl[:, :, round(Int64, t)]
     heatmap(frame_i, ratio=:equal, grid=false,
          xaxis=false, yaxis=false, xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0),
          colorbar_title="Voltage (mV)"
     )
end
gif(anim, "$(save_loc)/S6 ECl -55 Simulation.gif", fps=1000.0 / animate_dt)