using Revise
include("figure_setup.jl")
using RetinalChaos

import RetinalChaos.SDEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters
#Step 3: determine the timespan
tmin = 0.0
dt = 0.01
tmax = 120e3

#Step 4: set up the problem
parameters[I_app] = 1.0
parameters[g_Na] = 0.0
parameters[ρe] = 0.0
parameters[ρi] = 0.0
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
#Step 5: Solve the problem
@time sol = solve(probSDE, SOSRI()); #So far the best method is SOSRI

#First we want to show the role of Voltage gated sodium channels on spikes, bursts, and waves
#%% Step 2: Determine the number of trajectories and the parameter to adjust
reload_parameters()
parameters[I_app] = 1.0
parameters[g_Na] = 0.0
parameters[ρe] = 0.0
parameters[ρi] = 0.0
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
n_trajectories = 2
par = :g_Na #Adjust the noise
pmin = 0.0 #Determine the minimum parameter
pmax = 20.0 #Determine the maximum range of parameters
test_rng = LinRange(pmin, pmax, n_trajectories); #Determine the range of the parameters (specified above)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(probSDE, prob_func=prob_func); #Set up the problem
#Step 4: Run the simulation #Noise uses SOSRI(), 
@time sim = solve(ensemble_prob, SOSRA(), trajectories=n_trajectories, EnsembleThreads());

# Extract the data
cmap = [:red, :black]
#cmap = plt.get_cmap("plasma")
fig1, ax = plt.subplots(n_trajectories, 3, figsize=(10, 5))
for (sol_idx, sol_i) in enumerate(sim)
     tsSOL, dataSOL = timeseries_analysis(sol_i)
     vspikes = tsSOL["Spikes"][2]
     vbursts = tsSOL["Bursts"][2]
     println(vbursts)
     spike_start = vspikes[1, 1]
     spike_end = spike_start + 20.0
     tSpike = spike_start:dt:spike_end
     vSpike = map(t -> sol_i(t, idxs=2), tSpike)
     ax[sol_idx, 1].plot((tSpike .- tSpike[1]), vSpike, lw = 2.0, c=cmap[sol_idx])
     ax[sol_idx, 1].set_ylim(-85.0, 15.0)
     ax[sol_idx, 1].set_ylabel("Membrane \n Voltage (mV)", fontsize = 15.0)
     tBurst = vbursts[1, 1]-500:dt:vbursts[1, 1]+1500
     vBurst = map(t -> sol_i(t, idxs=2), tBurst)
     ax[sol_idx, 2].plot((tBurst .- tBurst[1])./1000, vBurst, lw = 2.0,  c=cmap[sol_idx])
     ax[sol_idx, 2].set_ylim(-85.0, 15.0)
     tIBI = sol.t
     vIBI = map(t -> sol_i(t, idxs=2), tIBI)
     ax[sol_idx, 3].plot((tIBI .- tIBI[1])./1000, vIBI, lw = 2.0, c=cmap[sol_idx])
     ax[sol_idx, 3].set_ylim(-85.0, 15.0)
end
ax[1, 1].xaxis.set_visible(false) #Turn off the bottom axis
ax[1, 1].legend(["gNa = 0.0 nS"], fontsize=12.0, handletextpad=0.5, loc = "lower center")
ax[2, 1].legend(["gNa = 2.0 nS"], fontsize=12.0, handletextpad=0.5, loc = "lower center")
ax[1, 1].spines["bottom"].set_visible(false)
ax[1, 2].xaxis.set_visible(false) #Turn off the bottom axis
ax[1, 2].spines["bottom"].set_visible(false)
ax[1, 3].xaxis.set_visible(false) #Turn off the bottom axis
ax[1, 3].spines["bottom"].set_visible(false)

ax[2, 1].set_xlabel("Time (ms)", fontsize = 15.0)
ax[2, 2].set_xlabel("Time (s)", fontsize = 15.0)
ax[2, 3].set_xlabel("Time (s)", fontsize = 15.0)

loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
print("[$(now())]: Saving the figure 2...")
fig1.savefig("$(loc)/S2 Sodium Conductances.jpg")
plt.close("all")
println(" Completed")