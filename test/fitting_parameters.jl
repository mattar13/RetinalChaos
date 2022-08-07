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
pars_dict[:C_m] = 13.6
#pars_dict[:I_app] = 10.0
pars_dict[:ρe] = 0.0
pars_dict[:ρi] = 0.0
pars_dict[:g_leak] = 1.718
p = pars_dict |> extract_dict
tspan = (0.0, 300e3)

#Questions to answer. 
#%% 1) How does the noiseless and noisy traces compare

#Noisy
probNOISE = SDEProblem(T_SDE, noise, u0, tspan, p)
@time solNOISE = solve(probNOISE, SOSRI(), save_idxs=1, saveat=dt, progress=true, progress_steps=1); #So far the best method is SOSRI
tsNOISE, dataNOISE = timeseries_analysis(solNOISE)

#No noise
probBASE = ODEProblem(T_ODE, u0, tspan, p)
@time solBASE = solve(probBASE, saveat=dt, save_idxs=1, progress=true, progress_steps=1); #So far the best method is SOSRI
tsBASE, dataBASE = timeseries_analysis(solBASE)

TimescaleLoss(data, solNOISE)

# Set the indexes for the spikes
#spike_t_BASE, spike_vt_BASE = RetinalChaos.extract_spike_trace(tsBASE, dataBASE, normalize=true)
spike_t_NOISE, spike_vt_NOISE = RetinalChaos.extract_spike_trace(tsNOISE, dataNOISE, normalize=true)

#burst_t_BASE, burst_vt_BASE = RetinalChaos.extract_burst_trace(tsBASE, dataBASE, normalize=true)
burst_t_NOISE, burst_vt_NOISE = RetinalChaos.extract_burst_trace(tsNOISE, dataNOISE, normalize=true)

#IBI_t_BASE, IBI_vt_BASE = RetinalChaos.extract_IBI_trace(tsBASE, dataBASE, normalize=true)
IBI_t_NOISE, IBI_vt_NOISE = RetinalChaos.extract_IBI_trace(tsNOISE, dataNOISE, normalize=true)

#%% calculate the loss on a single parameter set
TimescaleLoss(data, p; burst_weight=0.0, IBI_weight=0.0)
#%% Lets try to compute the gradient of each parameter 
fSPIKE(p) = TimescaleLoss(data, p; spike_weight=1.0, burst_weight=0.0, IBI_weight=0.0)
fBURST(p) = TimescaleLoss(data, p; spike_weight=0.0, burst_weight=1.0, IBI_weight=0.0)
fIBI(p) = TimescaleLoss(data, p; spike_weight=0.0, burst_weight=0.0, IBI_weight=1.0)
gSPIKE = x -> FD.gradient(fSPIKE, x)
gBURST = x -> FD.gradient(fBURST, x)
gIBI = x -> FD.gradient(fIBI, x)

MSE_grad_SPIKE = gSPIKE(p)
MSE_grad_BURST = gBURST(p)
MSE_grad_IBI = gIBI(p)

#Filter out only the highest changers
sort_spike = sortperm(MSE_grad_SPIKE)
sort_burst = sortperm(MSE_grad_BURST)
sort_IBI = sortperm(MSE_grad_IBI)

#%% Plot the comparison between the physiolgical trace
fig1, ax = plt.subplots(3, 3)
ax[1, 1].plot(spike_t_PHYS, spike_vt_PHYS, c=:black)
ax[2, 1].plot(spike_t_NOISE, spike_vt_NOISE, c=:blue)

ax[1, 2].plot(burst_t_PHYS, burst_vt_PHYS, c=:black)
ax[2, 2].plot(burst_t_NOISE, burst_vt_NOISE, c=:blue)

ax[1, 3].plot(IBI_t_PHYS, IBI_vt_PHYS, c=:black)
ax[2, 3].plot(IBI_t_NOISE, IBI_vt_NOISE, c=:blue)

ax[3, 1].barh(t_pars[sort_spike], MSE_grad_SPIKE[sort_spike])
ax[3, 2].barh(t_pars[sort_burst], MSE_grad_BURST[sort_burst])
ax[3, 3].barh(t_pars[sort_IBI], MSE_grad_IBI[sort_IBI])

#%% Why the heck is the delta parameter causing so much problems
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict #Initial conditions
pars_dict = read_JSON("params\\params.json")
#pars_dict[:I_app] = 5.0
pars_dict[:ρe] = 0.0
pars_dict[:ρi] = 0.0
p = pars_dict |> extract_dict #Parameters
tspan = (0.0, 300000.0) #Timespan
#prob = ODEProblem(T_ODE, u0, tspan, p) #ODE problem
prob = SDEProblem(T_SDE, noise, u0, tspan, p) #ODE problem

#Step 2: Determine the number of trajectories and the parameter to adjust
n_trajectories = 10
par_idx = p_find(:V2; list_p=t_pars); #Point to the index of the parameter
test_rng = LinRange(10.0, 20.0, n_trajectories); #Determine the range of the parameters (specified above)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #Set up the problem

#Step 4: Run the simulation #Noise uses SOSRI(), 
@time sim = solve(ensemble_prob, saveat=1.0, trajectories=n_trajectories, EnsembleThreads());

#%% Plot the solutions 
cmap = plt.get_cmap("plasma")
tsSIM, dataSIM = timeseries_analysis(sim[1])
fig, ax = plt.subplots(3, 3)
ax[1, 1].plot(spike_t_PHYS, spike_vt_PHYS, c=:black)
ax[1, 2].plot(burst_t_PHYS, burst_vt_PHYS, c=:black)
ax[1, 3].plot(IBI_t_PHYS, IBI_vt_PHYS, c=:black)

for (sol_idx, sol_i) in enumerate(sim)
     println(test_rng[sol_idx])
     tsSIM, dataSIM = timeseries_analysis(sol_i, idxs=1)
     if !isempty(tsSIM["Spikes"])
          spike_t_SIM, spike_vt_SIM = extract_spike_trace(tsSIM, dataSIM, idx=1, normalize=false)
          spike_nt_SIM = sol_i(spike_t_SIM, idxs=2) |> Array
          ax[2, 1].plot(spike_t_SIM .- spike_t_SIM[1], spike_vt_SIM, c=cmap(sol_idx / n_trajectories))
          ax[3, 1].plot(spike_nt_SIM, spike_vt_SIM, c=cmap(sol_idx / n_trajectories))
     end
     println(tsSIM["Bursts"])
     if !isempty(tsSIM["Bursts"])
          burst_t_SIM, burst_vt_SIM = extract_burst_trace(tsSIM, dataSIM, idx=1, normalize=false)
          burst_nt_SIM = sol_i(burst_t_SIM, idxs=2) |> Array

          IBI_t_SIM, IBI_vt_SIM = extract_IBI_trace(tsSIM, dataSIM, idx=1, normalize=false)
          IBI_nt_SIM = sol_i(IBI_t_SIM, idxs=2) |> Array

          ax[2, 2].plot(burst_t_SIM .- burst_t_SIM[1], burst_vt_SIM, c=cmap(sol_idx / n_trajectories))
          ax[3, 2].plot(burst_nt_SIM, burst_vt_SIM, c=cmap(sol_idx / n_trajectories))

          ax[2, 3].plot(IBI_t_SIM .- IBI_t_SIM[1], IBI_vt_SIM, c=cmap(sol_idx / n_trajectories))
          ax[3, 3].plot(IBI_nt_SIM, IBI_vt_SIM, c=cmap(sol_idx / n_trajectories))
     end
end
# creating ScalarMappable
sm = plt.cm.ScalarMappable(cmap=cmap)
cbar = fig.colorbar(sm, ticks = LinRange(0.0, 1.0, 4))
cbar.set_ticklabels(LinRange(test_rng[1], test_rng[end], 4))

#%%

#%%
idxs = map(x -> p_find(x), [:g_leak, :g_K, :g_Ca, :g_TREK, :C_m])
#Optimize only a few parameters
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

