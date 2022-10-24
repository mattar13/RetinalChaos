using Revise
include("figure_setup.jl")
using RetinalChaos

import RetinalChaos.SDEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters
#Step 3: determine the timespan
tmin = 0.0
dt = 0.01
tmax = 60e3

#Step 4: set up the problem
parameters[I_app] = 1.0
parameters[g_Na] = 0.0
parameters[ρe] = 0.0
parameters[ρi] = 0.0
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
#Step 5: Solve the problem
@time sol = solve(probSDE, SOSRI()); #So far the best method is SOSRI
tSOL = tmin:dt:tmax
ySOL = map(t -> sol(t, idxs = 2), tSOL)
p = plt.plot(tSOL, ySOL)
#%%

#First we want to show the role of Voltage gated sodium channels on spikes, bursts, and waves
#%% Step 2: Determine the number of trajectories and the parameter to adjust
reload_parameters()
parameters[I_app] = 1.0
parameters[g_Na] = 0.0
parameters[ρe] = 0.0
parameters[ρi] = 0.0
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
n_trajectories = 5
par = :g_Na #Adjust the noise
pmin = 0.0 #Determine the minimum parameter
pmax = 55.0 #Determine the maximum range of parameters
test_rng = LinRange(pmin, pmax, n_trajectories); #Determine the range of the parameters (specified above)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(probSDE, prob_func=prob_func); #Set up the problem
#Step 4: Run the simulation #Noise uses SOSRI(), 
@time sim = solve(ensemble_prob, SOSRA(), trajectories=n_trajectories, EnsembleThreads());

# Extract the data
cmap = plt.get_cmap("plasma")
fig1, ax = plt.subplots(3, 2, figsize=(10, 10))
for (sol_idx, sol_i) in enumerate(sim)
     tsSOL, dataSOL = timeseries_analysis(sol_i)
     vspikes = tsSOL["Spikes"][2]
     vbursts = tsSOL["Bursts"][2]
     println(vbursts)
     spike_start = vspikes[1, 1]
     spike_end = spike_start + 20.0
     tSpike = spike_start:dt:spike_end
     vSpike = map(t -> sol_i(t, idxs=2), tSpike)
     ax[1, 1].plot(tSpike .- tSpike[1], vSpike .+ (sol_idx * 10), c=cmap(sol_idx / n_trajectories))

     tBurst = vbursts[1, 1]-500:dt:vbursts[1, 1]+1500
     vBurst = map(t -> sol_i(t, idxs=2), tBurst)
     ax[2, 1].plot(tBurst .- tBurst[1], vBurst .+ (sol_idx * 50), c=cmap(sol_idx / n_trajectories))

     tIBI = vbursts[1, 1]-500:dt:vbursts[1, 1]+100e3
     vIBI = map(t -> sol_i(t, idxs=2), tIBI)
     ax[3, 1].plot(tIBI .- tIBI[1], vIBI .+ (sol_idx * 50), c=cmap(sol_idx / n_trajectories))
end
fig1