using Revise
include("figure_setup.jl")
using RetinalChaos

import RetinalChaos.SDEModel #import the ODEModel
import RetinalChaos.ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters

#%% How does each of the parameters affect the model statistics
using DataFrames, XLSX, Query
file = "$(loc)/change_data.xlsx"

#%% either we run the analysis again
#=
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
n = 80
percent = [ones(Int64(n / 2)); ones(Int64(n / 2)) .+ 0.1]
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
               if idx < Int64(n / 2)
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
          avg_spike0 = sum(spike_dur0) / length(spike_dur0) * 100
          avg_spike1 = sum(spike_dur1) / length(spike_dur1) * 100

          avg_burst0 = sum(burst_dur0) / length(burst_dur0) * 100
          avg_burst1 = sum(burst_dur1) / length(burst_dur1) * 100

          avg_IBI0 = sum(IBI_dur0) / length(IBI_dur0) * 100
          avg_IBI1 = sum(IBI_dur1) / length(IBI_dur1) * 100

          df[idx, :spike_percent] = (avg_spike0 / avg_spike1) - 1.0
          df[idx, :burst_percent] = (avg_burst0 / avg_burst1) - 1.0
          df[idx, :IBI_percent] = (avg_IBI0 / avg_IBI1) .- 1.0
     end
end
df = df |> @orderby_descending(_.params) |> DataFrame

rm(file)
XLSX.writetable(file, "PercentChange" => df)
=#
#%% Or we open the file

df = DataFrame(XLSX.readtable(file, "PercentChange"))
df.spike_percent

#%% Run a simulation
reload_parameters()
tmax = 120e3
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
@time sol = solve(probSDE, SOSRI(), save_idxs=2); #So far the best method is SOSRI
sol.u
timestamps, data = RetinalChaos.timeseries_analysis(sol)
tSPIKE, vSPIKE = RetinalChaos.extract_spike_trace(timestamps, data, normalize=false)
tBURST, vBURST = RetinalChaos.extract_burst_trace(timestamps, data, offset=500.0, burst_dur=1500, normalize=false)
tIBI, vIBI = RetinalChaos.extract_IBI_trace(timestamps, data, offset=1000, IBI_dur=60e3, normalize=false)
# We should pick a single parameter to demonstrate a 10% increase
fig2, ax = plt.subplots(2, 3, figsize=(7.5, 7.5),
     gridspec_kw=Dict("height_ratios" => (0.25, 0.75))
)
plt.subplots_adjust(
     left=0.17,
     right=0.90,
     bottom=0.12,
     top=0.9,
     wspace=0.5,
     hspace=0.25
)
col1_ylabel = -0.6
ax[1, 1].plot(tSPIKE .- tSPIKE[1], vSPIKE, c=:black, lw=1.0)
ax[1, 1].set_ylim(-90.0, 0.0)
ax[1, 1].set_xlabel("Time (ms)")
ax[1, 1].set_ylabel("Membrane Voltage (mV)")
ax[1, 1].yaxis.set_label_coords(col1_ylabel, 0.5)
ax[1, 1].yaxis.set_major_locator(MultipleLocator(40.0))
ax[1, 1].yaxis.set_minor_locator(MultipleLocator(20.0))
ax[1, 1].xaxis.set_major_locator(MultipleLocator(12.0))
ax[1, 1].xaxis.set_minor_locator(MultipleLocator(6.0))
ax[1, 1].set_title("Spike")

ax[1, 2].plot((tBURST .- tBURST[1]) ./ 1000, vBURST, c=:black, lw=1.0)
ax[1, 2].set_ylim(-90.0, 0.0)
ax[1, 2].set_xlabel("Time (s)")
ax[1, 2].yaxis.set_major_locator(MultipleLocator(40.0))
ax[1, 2].yaxis.set_minor_locator(MultipleLocator(20.0))
ax[1, 2].xaxis.set_major_locator(MultipleLocator(1.0))
ax[1, 2].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[1, 2].yaxis.set_visible(false) #Turn off the bottom axis
ax[1, 2].spines["left"].set_visible(false)
ax[1, 2].set_title("Burst")

ax[1, 3].plot((tIBI .- tIBI[1]) ./ 1000, vIBI, c=:black, lw=1.0)
ax[1, 3].set_ylim(-90.0, 0.0)
ax[1, 3].set_xlabel("Time (s)")
ax[1, 3].yaxis.set_major_locator(MultipleLocator(40.0))
ax[1, 3].yaxis.set_minor_locator(MultipleLocator(20.0))
ax[1, 3].xaxis.set_major_locator(MultipleLocator(30.0))
ax[1, 3].xaxis.set_minor_locator(MultipleLocator(15.0))
ax[1, 3].yaxis.set_visible(false) #Turn off the bottom axis
ax[1, 3].spines["left"].set_visible(false)
ax[1, 3].set_title("Interburst Interval")

plot_params = df.params[findall(df.params .∉ Ref(["g_ACh", "g_GABA", "ρe", "ρi", "E_ACh", "E_Cl", "τ_ACh", "τ_GABA", "k_GABA", "k_ACh"]))]

dfSPIKE = df |> @orderby(_.spike_percent) |> @filter(_.params ∈ plot_params) |> DataFrame
nshow = 5
ax[2, 1].barh(dfSPIKE.params[1:nshow], (dfSPIKE.spike_percent[1:nshow]) .* 100, color=:red, height=0.5)
ax[2, 1].barh(dfSPIKE.params[end-(nshow-1):end], (dfSPIKE.spike_percent[end-(nshow-1):end]) .* 100, color=:blue, height=0.5)
ax[2, 1].vlines([0.0], ymin=-0.5, ymax=(nshow * 2.0)-0.5, color=:black, lw=2.0)
ax[2, 1].set_xlabel("Spike Duration \n Increase (%)")
ax[2, 1].set_ylabel("Parameter 10% Increase")
ax[2, 1].yaxis.set_label_coords(col1_ylabel, 0.5)

dfBURST = df |> @orderby(_.burst_percent) |> @filter(_.params ∈ plot_params) |> DataFrame
ax[2, 2].barh(dfBURST.params[1:nshow], (dfBURST.burst_percent[1:nshow]) .* 100, color=:red, height=0.5)
ax[2, 2].barh(dfBURST.params[end-(nshow-1):end], (dfBURST.burst_percent[end-(nshow-1):end]) .* 100, color=:blue, height=0.5)
ax[2, 2].vlines([0.0], ymin=-0.5, ymax=(nshow * 2.0)-0.5, color=:black, lw=2.0)
ax[2, 2].set_xlabel("Burst Duration \n Increase (%)")

dfIBI = df |> @orderby(_.IBI_percent) |> @filter(_.params ∈ plot_params) |> DataFrame
ax[2, 3].barh(dfIBI.params[1:nshow], (dfIBI.IBI_percent[1:nshow]) .* 100, color=:red, height=0.5)
ax[2, 3].barh(dfIBI.params[end-(nshow-1):end], (dfIBI.IBI_percent[end-(nshow-1):end]) .* 100, color=:blue, height=0.5)
ax[2, 3].vlines([0.0], ymin=-0.5, ymax=(nshow * 2.0)-0.5, color=:black, lw=2.0)
ax[2, 3].set_xlabel("IBI \n Increase (%)")

ax[2, 3].annotate("A", (0.01, 0.94), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
ax[2, 3].annotate("B", (0.01, 0.60), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
#%%
print("[$(now())]: Saving the figure 2...")
fig2.savefig("$(loc)/Figure4_Altering_Parameters.jpg")
plt.close("all")
println(" Completed")


#%% 
reload_parameters()
n = 20
par = :V2
rng = LinRange(20.0, 30.0, n)
tmax = 2e3
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par, rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(probSDE, prob_func=prob_func) #Set up the problem
@time sim = solve(ensemble_prob, SOSRI(), trajectories=n, EnsembleThreads())

fig2, ax = plt.subplots(1)
for sol in sim
     tSERIES = sol.t
     xSERIES = map(t -> sol(t)[2], tSERIES)
     ax.plot(tSERIES, xSERIES)
end
fig2