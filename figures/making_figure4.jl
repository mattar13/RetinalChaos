
#%% How does each of the parameters affect the model statistics
using DataFrames, XLSX, Query

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

rm("$(loc)/change_data.xlsx")
XLSX.writetable("$(loc)/change_data.xlsx", "PercentChange" => df)

#%% Run a simulation
reload_parameters()
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
@time sol = solve(probSDE, SOSRI(), save_idxs=2); #So far the best method is SOSRI
sol
timestamps, data = RetinalChaos.timeseries_analysis(sol)
tSPIKE, vSPIKE = RetinalChaos.extract_spike_trace(timestamps, data, normalize=false)
tBURST, vBURST = RetinalChaos.extract_burst_trace(timestamps, data, offset = 500.0, burst_dur=1500, normalize=false)
tIBI, vIBI = RetinalChaos.extract_IBI_trace(timestamps, data, offset = 1000, IBI_dur = 60e3, normalize=false)
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

ax[1, 2].plot((tBURST .- tBURST[1]) ./ 1000, vBURST, c=:black, lw=1.0)
ax[1, 2].set_ylim(-90.0, 0.0)
ax[1, 2].set_xlabel("Time (s)")
ax[1, 2].yaxis.set_major_locator(MultipleLocator(40.0))
ax[1, 2].yaxis.set_minor_locator(MultipleLocator(20.0))
ax[1, 2].xaxis.set_major_locator(MultipleLocator(1.0))
ax[1, 2].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[1, 2].yaxis.set_visible(false) #Turn off the bottom axis
ax[1, 2].spines["left"].set_visible(false)

ax[1, 3].plot((tIBI .- tIBI[1]) ./ 1000, vIBI, c=:black, lw=1.0)
ax[1, 3].set_ylim(-90.0, 0.0)
ax[1, 3].set_xlabel("Time (s)")
ax[1, 3].yaxis.set_major_locator(MultipleLocator(40.0))
ax[1, 3].yaxis.set_minor_locator(MultipleLocator(20.0))
ax[1, 3].xaxis.set_major_locator(MultipleLocator(30.0))
ax[1, 3].xaxis.set_minor_locator(MultipleLocator(15.0))
ax[1, 3].yaxis.set_visible(false) #Turn off the bottom axis
ax[1, 3].spines["left"].set_visible(false)

dfSPIKE = df |> @orderby(_.spike_percent) |> @filter(!isnan(_.spike_percent)) |> DataFrame
ax[2, 1].barh(dfSPIKE.params[1:10], (dfSPIKE.spike_percent[1:10]) .* 100, color=:red, height=0.5)
ax[2, 1].barh(dfSPIKE.params[end-9:end], (dfSPIKE.spike_percent[end-9:end]) .* 100, color=:blue, height=0.5)
ax[2, 1].vlines([0.0], ymin=-1.0, ymax=20.0, color=:black, lw=2.0)
ax[2, 1].set_xlabel("Spike Duration \n Increase (%)")
ax[2, 1].set_ylabel("Parameter 10% Increase")
ax[2, 1].yaxis.set_label_coords(col1_ylabel, 0.5)

dfBURST = df |> @orderby(_.burst_percent) |> @filter(!isnan(_.burst_percent)) |> DataFrame
ax[2, 2].barh(dfBURST.params[1:10], (dfBURST.burst_percent[1:10]) .* 100, color=:red, height=0.5)
ax[2, 2].barh(dfBURST.params[end-9:end], (dfBURST.burst_percent[end-9:end]) .* 100, color=:blue, height=0.5)
ax[2, 2].vlines([0.0], ymin=-1.0, ymax=20.0, color=:black, lw=2.0)
ax[2, 2].set_xlabel("Burst Duration \n Increase (%)")

dfIBI = df |> @orderby(_.IBI_percent) |> @filter(!isnan(_.IBI_percent)) |> DataFrame
ax[2, 3].barh(dfIBI.params[1:10], (dfIBI.IBI_percent[1:10]) .* 100, color=:red, height=0.5)
ax[2, 3].barh(dfIBI.params[end-9:end], (dfIBI.IBI_percent[end-9:end]) .* 100, color=:blue, height=0.5)
ax[2, 3].vlines([0.0], ymin=-1.0, ymax=20.0, color=:black, lw=2.0)
ax[2, 3].set_xlabel("IBI \n Increase (%)")

#%%
print("[$(now())]: Saving the figure 2...")
fig2.savefig("$(loc)/Figure4_Altering_Parameters.jpg")
plt.close("all")
println(" Completed")