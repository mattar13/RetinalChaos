using Revise
using RetinalChaos
using DataFrames
include("figure_setup.jl")

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using ProgressMeter

using Plots
# Set some initial parameters for this script
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling"
nx = ny = 64
nt = 120001
animate_dt = 60.0

#%% Model 1: Regular Baseline model 
wave_path = "$(data_root)\\wave_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
@time res = run_model(p_dict, u_dict, wave_path; 
     run_simulation = false, burst_wave = false, verbose = true, 
     alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7
)

#%% Model 2: Blocked Neurotransmission 
isolated_path = "$(data_root)/isolated_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:g_ACh] = 0.0 # Block all Acetylcholine receptors
p_dict[:g_GABA] = 0.0 #Block all GABA receptors
@time res = run_model(p_dict, u_dict, isolated_path; 
     run_simulation = false, burst_wave = false, verbose = true, 
     alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7
)

#%% Model 3: No GABA
noGABA_path = "$(data_root)/no_GABA_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:g_GABA] = 0.0 #Block all GABA receptors
@time res = run_model(p_dict, u_dict, noGABA_path; 
     run_simulation = true, burst_wave = true, verbose = true, 
     alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7
)

#%% Model 4 ECl Differential
ec_path = "$(data_root)\\ECl55_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:E_Cl] = -55.0 #Block all GABA receptors
@time res = run_model(p_dict, u_dict, ec_path; 
     run_simulation = true, burst_wave = true, verbose = true, 
     alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7
)

#%% Model 5 ECl more hyperpolarizing Differential
ecH_path = "$(data_root)\\ECl75_model"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:E_Cl] = -75.0 #Block all GABA receptors
@time res = run_model(p_dict, u_dict, ecH_path; 
     run_simulation = true, burst_wave = true, verbose = true, 
     alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7
)

#Open up data for the Diffusion experiments==================================================================================================================#

#%% What are the effects of changing the synaptic release? 
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\FastEdiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:De] = 0.01 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#Fast inhibitory ==================================================================================================================#
#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\FastIDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:Di] = 0.01 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\FastI_EDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:De] = 0.01 #Make synaptic transmission faster
p_dict[:Di] = 0.01 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% What are the effects of changing the synaptic release? 
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\SlowEdiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:De] = 0.001 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\SlowIDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:Di] = 0.001 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)

#%% Change the inhibitory diffusion properties
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\figure_data\SlowI_EDiffusion"
u_dict = read_JSON("params/conds.json")
p_dict = read_JSON("params/params.json")
p_dict[:Di] = 0.001 #Make synaptic transmission faster
p_dict[:De] = 0.001 #Make synaptic transmission faster
run_model(u_dict, p_dict, loc; alg = SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7)