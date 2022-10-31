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
     ax[sol_idx, 1].plot((tSpike .- tSpike[1]), vSpike, lw=2.0, c=cmap[sol_idx])
     ax[sol_idx, 1].set_ylim(-85.0, 15.0)
     ax[sol_idx, 1].set_ylabel("Membrane \n Voltage (mV)", fontsize=15.0)
     tBurst = vbursts[1, 1]-500:dt:vbursts[1, 1]+1500
     vBurst = map(t -> sol_i(t, idxs=2), tBurst)
     ax[sol_idx, 2].plot((tBurst .- tBurst[1]) ./ 1000, vBurst, lw=2.0, c=cmap[sol_idx])
     ax[sol_idx, 2].set_ylim(-85.0, 15.0)
     tIBI = sol.t
     vIBI = map(t -> sol_i(t, idxs=2), tIBI)
     ax[sol_idx, 3].plot((tIBI .- tIBI[1]) ./ 1000, vIBI, lw=2.0, c=cmap[sol_idx])
     ax[sol_idx, 3].set_ylim(-85.0, 15.0)
end
ax[1, 1].xaxis.set_visible(false) #Turn off the bottom axis
ax[1, 1].legend(["gNa = 0.0 nS"], fontsize=12.0, handletextpad=0.5, loc="lower center")
ax[2, 1].legend(["gNa = 2.0 nS"], fontsize=12.0, handletextpad=0.5, loc="lower center")
ax[1, 1].spines["bottom"].set_visible(false)
ax[1, 2].xaxis.set_visible(false) #Turn off the bottom axis
ax[1, 2].spines["bottom"].set_visible(false)
ax[1, 3].xaxis.set_visible(false) #Turn off the bottom axis
ax[1, 3].spines["bottom"].set_visible(false)

ax[2, 1].set_xlabel("Time (ms)", fontsize=15.0)
ax[2, 2].set_xlabel("Time (s)", fontsize=15.0)
ax[2, 3].set_xlabel("Time (s)", fontsize=15.0)

loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
print("[$(now())]: Saving the figure 2...")
fig1.savefig("$(loc)/S2 Sodium Conductances.jpg")
plt.close("all")
println(" Completed")

#%% How does each of the parameters affect the model statistics
using DataFrames, XLSX, Query

probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
n = 80
percent = [ones(Int64(n/2)); ones(Int64(n/2)).+0.1]
n_trajectories = length(percent)
df = DataFrame(params="", initial=ones(length(parameters)), increase=1.0, spike_percent=1.0, burst_percent=1.0, IBI_percent=1.0)

for (idx, entry) in enumerate(parameters)
     reload_parameters()
     par = Symbol(entry[1])
     println(par)
     df[idx, :params] = String(par)
     df[idx, :initial] = entry[2]
     df[idx, :increase] = entry[2] * 1.1
     test_rng = entry[2] .* percent #Determine the range of the parameters (specified above)
     if par ∈ [:Di, :De, :σ, :I_app]
          df[idx, :spike_percent] = 0.0
          df[idx, :burst_percent] = 0.0
          df[idx, :IBI_percent] = 0.0
     else
          prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par, test_rng) #Set up the problem function to indicate that the voltage will be altered
          ensemble_prob = EnsembleProblem(probSDE, prob_func=prob_func) #Set up the problem
          @time sim = solve(ensemble_prob, SOSRA(), trajectories=n_trajectories, EnsembleThreads())
          spike_dur0 = Float64[]
          burst_dur0 = Float64[]
          IBI_dur0 = Float64[]

          spike_dur1 = Float64[]
          burst_dur1 = Float64[]
          IBI_dur1 = Float64[]

          for (idx, sol) in enumerate(sim)
               timestamps, data = timeseries_analysis(sol)
               if idx < Int64(n/2)
                    if !isnan(data["SpikeDurAvg"])
                         push!(spike_dur0, data["SpikeDurAvg"])
                    end
                    if !isnan(data["BurstDurAvg"])
                         push!(burst_dur0, data["BurstDurAvg"])
                    end
                    if !isnan(data["IBIAvg"])
                         push!(IBI_dur0, data["IBIAvg"])
                    end
               else
                    if !isnan(data["SpikeDurAvg"])
                         push!(spike_dur1, data["SpikeDurAvg"])
                    end
                    if !isnan(data["BurstDurAvg"])
                         push!(burst_dur1, data["BurstDurAvg"])
                    end
                    if !isnan(data["IBIAvg"])
                         push!(IBI_dur1, data["IBIAvg"])
                    end
               end
               #println(keys(data))
          end
          avg_spike0 = sum(spike_dur0)/length(spike_dur0)
          avg_spike1 = sum(spike_dur1)/length(spike_dur1)

          avg_burst0 = sum(burst_dur0)/length(burst_dur0)
          avg_burst1 = sum(burst_dur1)/length(burst_dur1)

          avg_IBI0 = sum(IBI_dur0)/length(IBI_dur0)
          avg_IBI1 = sum(IBI_dur1)/length(IBI_dur1)

          df[idx, :spike_percent] = (avg_spike0 / avg_spike1) - 1.0
          df[idx, :burst_percent] = (avg_burst0 / avg_burst1) - 1.0
          df[idx, :IBI_percent] = (avg_IBI0 / avg_IBI1) .- 1.0
     end
end
df = df |> @orderby_descending(_.params) |> DataFrame

rm("$(loc)/change_data.xlsx")
XLSX.writetable("$(loc)/change_data.xlsx", "PercentChange" => df)

fig2, ax = plt.subplots(3, figsize=(7.5, 7.5))
plt.subplots_adjust(
     left=0.15,
     right=0.95,
     bottom=0.12,
     top=0.9,
     wspace=0.0,
     hspace=0.7
)

dfSPIKE = df |> @orderby_descending(_.spike_percent) |> @filter(!isnan(_.spike_percent)) |> DataFrame
ax[1].bar(dfSPIKE.params[1:10], (dfSPIKE.spike_percent[1:10]) .* 100, color=:blue)
ax[1].bar(dfSPIKE.params[end-9:end], (dfSPIKE.spike_percent[end-9:end]) .* 100, color=:red)
#ax[1].set_ylim(-25.0, 100.0)
ax[1].hlines([0.0], xmin=-0.5, xmax=9.5, color=:black)
ax[1].tick_params(axis="x", rotation=90)
ax[1].set_title("Parameter 10% Increase")
ax[1].set_ylabel("Spike \n Percent Increase")

dfBURST = df |> @orderby_descending(_.burst_percent) |> @filter(!isnan(_.burst_percent)) |> DataFrame
ax[2].bar(dfBURST.params[1:10], (dfBURST.burst_percent[1:10]) .* 100, color=:blue)
ax[2].bar(dfBURST.params[end-9:end], (dfBURST.burst_percent[end-9:end]) .* 100, color=:red)
#ax[2].set_ylim(-75.0, 50.0)
ax[2].hlines([0.0], xmin=-0.5, xmax=9.5, color=:black)
ax[2].tick_params(axis="x", rotation=90)
ax[2].set_ylabel("Burst \n Percent Increase")

dfIBI = df |> @orderby_descending(_.IBI_percent) |> @filter(!isnan(_.IBI_percent)) |> DataFrame
ax[3].bar(dfIBI.params[1:10], (dfIBI.IBI_percent[1:10]) .* 100, color=:blue)
ax[3].bar(dfIBI.params[end-9:end], (dfIBI.IBI_percent[end-9:end]) .* 100, color=:red)
#ax[3].set_ylim(-75.0, 150.0)
ax[3].hlines([0.0], xmin=-0.5, xmax=9.5, color=:black)
ax[3].tick_params(axis="x", rotation=90)
ax[3].set_ylabel("IBI \n Percent Increase")

print("[$(now())]: Saving the figure 2...")
fig2.savefig("$(loc)/S7 Sodium Conductances.jpg")
plt.close("all")
println(" Completed")
