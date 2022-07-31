using Revise
using RetinalChaos
using ePhys
using Statistics, StatsBase
import RetinalChaos: MeanSquaredError
using Optim
include("../figures/figure_setup.jl") #This will load all of the data we want
include("../figures/opening_data.jl") #This will load all of the data we want

#%% Set up the parameters we will use throught the fitting
dt = 0.5
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 10.0
pars_dict[:ρe] = 0.0
pars_dict[:ρi] = 0.0
p = pars_dict |> extract_dict
tspan = (0.0, 300e3)

#Questions to answer. 
#%% 1) How does the noiseless and noisy traces compare

#Noisy
probNOISE = SDEProblem(T_SDE, noise, u0, tspan, p)
@time solNOISE = solve(probNOISE, SOSRI(), saveat=dt, progress=true, progress_steps=1); #So far the best method is SOSRI
tsNOISE, dataNOISE = timeseries_analysis(solNOISE)

#No noise
probBASE = ODEProblem(T_ODE, u0, tspan, p)
@time solBASE = solve(probBASE, saveat=dt, progress=true, progress_steps=1); #So far the best method is SOSRI
tsBASE, dataBASE = timeseries_analysis(solBASE)
TimescaleLoss(solBASE, solNOISE, dt=dt)

# Set the indexes for the spikes
spike_t_BASE, spike_vt_BASE = RetinalChaos.extract_spike_trace(tsBASE, dataBASE, normalize = true)
spike_t_NOISE, spike_vt_NOISE = RetinalChaos.extract_spike_trace(tsNOISE, dataNOISE, normalize = true)
spike_MSE = MeanSquaredError(spike_vt_BASE, spike_vt_NOISE)

burst_t_BASE, burst_vt_BASE = RetinalChaos.extract_burst_trace(tsBASE, dataBASE, normalize = true)
burst_t_NOISE, burst_vt_NOISE = RetinalChaos.extract_burst_trace(tsNOISE, dataNOISE, normalize = true)
burst_MSE = MeanSquaredError(burst_vt_BASE, burst_vt_NOISE)

IBI_t_BASE, IBI_vt_BASE = RetinalChaos.extract_IBI_trace(tsBASE, dataBASE, normalize = true)
IBI_t_NOISE, IBI_vt_NOISE = RetinalChaos.extract_IBI_trace(tsNOISE, dataNOISE, normalize = true)
IBI_MSE = MeanSquaredError(IBI_vt_BASE, IBI_vt_NOISE)

fig1, ax = plt.subplots(3, 3)
ax[1, 1].plot(spike_t_BASE, spike_vt_BASE)
ax[2, 1].plot(spike_t_NOISE, spike_vt_NOISE)

ax[1, 2].plot(burst_t_BASE, burst_vt_BASE)
ax[2, 2].plot(burst_t_NOISE, burst_vt_NOISE)

ax[1, 3].plot(IBI_t_BASE, IBI_vt_BASE)
ax[2, 3].plot(IBI_t_NOISE, IBI_vt_NOISE)

ax[3, 1].bar(["Spike", "Burst", "IBI"], [spike_MSE, burst_MSE, IBI_MSE])

#%% 2) Compare the simulation with the real data
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
#Eventually we will open multiple analysis files
data = readABF(target_file, channels=["Vm_prime4"], stimulus_name=nothing, time_unit=:ms)
downsample!(data, 1 / dt) #It is necessary to downsample to get understandable results
tsPHYS, dataPHYS = timeseries_analysis(data) #Move timeseries analysis to ePhys
TimescaleLoss(solNOISE, data)
#%% 
spike_t_PHYS, spike_vt_PHYS = RetinalChaos.extract_spike_trace(tsPHYS, dataPHYS, idx=3, dt=0.5)
burst_t_PHYS, burst_vt_PHYS = RetinalChaos.extract_burst_trace(tsPHYS, dataPHYS, idx=2, dt=0.5)
IBI_t_PHYS, IBI_vt_PHYS = RetinalChaos.extract_IBI_trace(tsPHYS, dataPHYS, idx=2, dt=0.5)

fig2, ax = plt.subplots(3, 3)
ax[1, 1].plot(spike_t_PHYS, tf)
ax[2, 1].plot(spike_t_NOISE, spike_vt_NOISE./maximum(spike_vt_NOISE))

ax[1, 2].plot(burst_t_PHYS, burst_vt_PHYS)
ax[2, 2].plot(burst_t_NOISE, burst_vt_NOISE)

ax[1, 3].plot(IBI_t_PHYS, IBI_vt_PHYS)
ax[2, 3].plot(IBI_t_NOISE, IBI_vt_NOISE)

TimescaleLoss(data, p; burst_weight=0.0, IBI_weight=0.0)
#%% Lets try to compute the gradient
fMIN(p) = TimescaleLoss(data, p)
gMIN = x -> FD.gradient(fMIN, x)
gradMSE = gMIN(p)
sort_idxs = sortperm(gradMSE)
fig, ax = plt.subplots(1)
ax.bar(t_pars[sort_idxs], gradMSE[sort_idxs])
t_pars
#%%
#Optimize only a few parameters
idxs = map(x -> p_find(x), [:g_leak, :g_K, :g_Ca, :g_TREK, :C_m])
function fSPIKE(data, pSECTION; test_parameters=[:g_leak, :g_K, :g_Ca, :g_TREK, :C_m])
     pars_dict = read_JSON("params\\params.json")
     p = pars_dict |> extract_dict
     idxs = map(x -> p_find(x), test_parameters)
     p[idxs] .= pSECTION
     TimescaleLoss(data, p; burst_weight=0.0, IBI_weight=0.0)
end

p_test = p[idxs]
fSPIKE(data, p_test)

results = optimize(fMIN, p)
pOPT = Optim.minimizer(results)


#%%
fig, ax = plt.subplots(1)
ax.plot(losses)
argmin(losses)
#%% Try plotting the improved function
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
pOPT = param_set[:, argmin(losses)]
tspan = (0.0, 300e3)
#Noisy
probOPT = SDEProblem(T_SDE, noise, u0, tspan, pOPT)
@time solOPT = solve(probOPT, SOSRI(), saveat=dt, progress=true, progress_steps=1); #So far the best method is SOSRI
plot(solOPT)
tsOPT, dataOPT = timeseries_analysis(solOPT)
spike_t_OPT, spike_vt_OPT = RetinalChaos.extract_spike_trace(tsOPT, dataOPT)
burst_t_OPT, burst_vt_OPT = RetinalChaos.extract_burst_trace(tsOPT, dataOPT)
IBI_t_OPT, IBI_vt_OPT = RetinalChaos.extract_IBI_trace(tsOPT, dataOPT)

fig, ax = plt.subplots(3, 3)
ax[1, 1].plot(spike_t_PHYS, spike_vt_PHYS)
ax[2, 1].plot(spike_t_OPT, spike_vt_OPT)

ax[1, 2].plot(burst_t_PHYS, burst_vt_PHYS)
ax[2, 2].plot(burst_t_OPT, burst_vt_OPT)

ax[1, 3].plot(IBI_t_PHYS, IBI_vt_PHYS)
ax[2, 3].plot(IBI_t_OPT, IBI_vt_OPT)

