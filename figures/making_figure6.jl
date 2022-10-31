using Revise
include("figure_setup.jl")
using Dates
using RetinalChaos
using StatsBase, Statistics
using ePhys
import ePhys: dwt_filter
include("opening_data.jl")

#%%
print("[$(now())]: Plotting... ")
width_inches = 7.5
height_inches = 7.5
fig4 = plt.figure("Physiology Data", figsize=(width_inches, height_inches))

col1_ylabel = -0.21
col2_ylabel = -0.11

gs = fig4.add_gridspec(5, 3,
     #width_ratios=(0.50, 0.50),
     #height_ratios=(0.16, 0.16, 0.15, 0.15, 0.15, 0.15),
     right=0.95, left=0.15,
     top=0.95, bottom=0.08,
     wspace=0.20, hspace=0.5)

#axA = fig4.add_subplot(py"""$(gs)[0, :]""")
axA1 = fig4.add_subplot(gs[1, 1])
ylim(-90.0, 10.0)
ylabel("Whole Cell \n Vm (mV)")
axA1.set_title("Spikes")
axA1.plot(spike_t_PHYS ./ 1000, spike_vt_PHYS, c=:black, lw=lw_standard)
axA1.xaxis.set_visible(false) #Turn off the bottom axis
axA1.spines["bottom"].set_visible(false)
axA1.yaxis.set_label_coords(col1_ylabel, 0.5)
axA1.yaxis.set_major_locator(MultipleLocator(20.0))
axA1.yaxis.set_minor_locator(MultipleLocator(10.0))

axA2 = fig4.add_subplot(gs[1, 2])
ylim(-90.0, 10.0)
axA2.set_title("Bursts")
axA2.plot(burst_t_PHYS ./ 1000, burst_vt_PHYS, c=:black, lw=lw_standard)
axA2.xaxis.set_visible(false) #Turn off the bottom axis
axA2.spines["bottom"].set_visible(false)
axA2.yaxis.set_major_locator(MultipleLocator(40.0))
axA2.yaxis.set_minor_locator(MultipleLocator(20.0))

axA3 = fig4.add_subplot(gs[1, 3])
ylim(-90.0, 10.0)
axA3.set_title("Interburst Interval")
axA3.plot(IBI_t_PHYS ./ 1000, IBI_vt_PHYS, c=:black, lw=lw_standard)
axA3.xaxis.set_visible(false)
axA3.spines["bottom"].set_visible(false)
axA3.yaxis.set_major_locator(MultipleLocator(40.0))
axA3.yaxis.set_minor_locator(MultipleLocator(20.0))

# Figure Part B
#axB = fig4.add_subplot(py"""$(gs)[1, :]""")
axB1 = fig4.add_subplot(gs[2, 1])
ylim(-90.0, 10.0)
ylabel("Isolated \n Vt (mV)")
axB1.plot(spike_t_ISO ./ 1000, spike_vt_ISO, c=:red, lw=lw_standard)
axB1.xaxis.set_visible(false) #Turn off the bottom axis
axB1.spines["bottom"].set_visible(false)
axB1.yaxis.set_label_coords(col1_ylabel, 0.5)
axB1.yaxis.set_major_locator(MultipleLocator(20.0))
axB1.yaxis.set_minor_locator(MultipleLocator(10.0))

axB2 = fig4.add_subplot(gs[2, 2])
ylim(-90.0, 10.0)
axB2.plot(burst_t_ISO ./ 1000, burst_vt_ISO, c=:red, lw=lw_standard)
axB2.xaxis.set_visible(false) #Turn off the bottom axis
axB2.spines["bottom"].set_visible(false)
axB2.yaxis.set_major_locator(MultipleLocator(40.0))
axB2.yaxis.set_minor_locator(MultipleLocator(20.0))

axB3 = fig4.add_subplot(gs[2, 3])
ylim(-90.0, 10.0)
axB3.plot(IBI_t_ISO ./ 1000, IBI_vt_ISO, c=:red, lw=lw_standard)
axB3.xaxis.set_visible(false) #Turn off the bottom axis
axB3.spines["bottom"].set_visible(false)
axB3.yaxis.set_major_locator(MultipleLocator(40.0))
axB3.yaxis.set_minor_locator(MultipleLocator(20.0))

# Figure part C Normal Wave model
axC1 = fig4.add_subplot(gs[3, 1])
ylim(-90.0, 10.0)
ylabel("ECl -55 \n Vt (mV)")
axC1.plot(spike_t_EC ./ 1000, spike_vt_EC, c=:blue, lw=lw_standard)
axC1.xaxis.set_visible(false) #Turn off the bottom axis
axC1.spines["bottom"].set_visible(false)
axC1.yaxis.set_label_coords(col1_ylabel, 0.5)
axC1.yaxis.set_major_locator(MultipleLocator(20.0))
axC1.yaxis.set_minor_locator(MultipleLocator(10.0))

axC2 = fig4.add_subplot(gs[3, 2])
ylim(-90.0, 10.0)
axC2.plot(burst_t_EC, burst_vt_EC, c=:blue, lw=lw_standard)
axC2.xaxis.set_visible(false) #Turn off the bottom axis
axC2.spines["bottom"].set_visible(false)
axC2.yaxis.set_major_locator(MultipleLocator(40.0))
axC2.yaxis.set_minor_locator(MultipleLocator(20.0))

#Need to run this one again
axC3 = fig4.add_subplot(gs[3, 3])
ylim(-90.0, 10.0)
axC3.plot(IBI_t_EC ./ 1000, IBI_vt_EC, c=:blue, lw=lw_standard)
axC3.xaxis.set_visible(false) #Turn off the bottom axis
axC3.spines["bottom"].set_visible(false)
axC3.yaxis.set_major_locator(MultipleLocator(40.0))
axC3.yaxis.set_minor_locator(MultipleLocator(20.0))

#axC = fig4.add_subplot(py"""$(gs)[2, :]""")
axD1 = fig4.add_subplot(gs[4, 1])
ylim(-90.0, 10.0)
ylabel("ECl -65 \n Vt (mV)")
xlabel("Time (s)")
axD1.plot((spike_t_WAVE.-spike_t_WAVE[1]) ./ 1000, spike_vt_WAVE, c=:green, lw=lw_standard)
axD1.yaxis.set_label_coords(col1_ylabel, 0.5)
axD1.yaxis.set_major_locator(MultipleLocator(20.0))
axD1.yaxis.set_minor_locator(MultipleLocator(10.0))

axD2 = fig4.add_subplot(gs[4, 2])
ylim(-90.0, 10.0)
xlabel("Time (s)")
axD2.plot((burst_t_WAVE.-burst_t_WAVE[1]) ./ 1000, burst_vt_WAVE, c=:green, lw=lw_standard)
axD2.yaxis.set_major_locator(MultipleLocator(40.0))
axD2.yaxis.set_minor_locator(MultipleLocator(20.0))

axD3 = fig4.add_subplot(gs[4, 3])
ylim(-90.0, 10.0)
xlabel("Time (s)")
axD3.plot((IBI_t_WAVE.-IBI_t_WAVE[1]) ./ 1000, IBI_vt_WAVE, c=:green, lw=lw_standard)
axD3.yaxis.set_major_locator(MultipleLocator(40.0))
axD3.yaxis.set_minor_locator(MultipleLocator(20.0))
#Make Bar Graphs showing the MSE loss between each RandomNumbers

jitter_dataPHYS_spike = 0.5 .* rand(length(phys_spike_durs))
jitter_dataISO_spike = 0.5 .* rand(length(dataISO["SpikeDurs"])) .+ 1
jitter_dataEC_spike = 0.5 .* rand(length(dataEC["SpikeDurs"])) .+ 2
jitter_dataWAVE_spike = 0.5 .* rand(length(dataWAVE["SpikeDurs"])) .+ 3

labels = ["PHYS", "ISO", "ECl", "WAV"]
pos = [0.25, 1.25, 2.25, 3.25]
axE1 = fig4.add_subplot(gs[5, 1])
#ylim(0.0, 100.0)
#axE1.plot(sdur_weights, sdur_edges, c = :black)
axE1.scatter(jitter_dataPHYS_spike, phys_spike_durs./1000, s=5.0, c=:black)
axE1.scatter(jitter_dataISO_spike, dataISO["SpikeDurs"]./1000, s=5.0, c=:red)
axE1.scatter(jitter_dataEC_spike, dataEC["SpikeDurs"]./1000, s=5.0, c=:blue)
axE1.scatter(jitter_dataWAVE_spike, dataWAVE["SpikeDurs"]./1000, s=5.0, c=:green)
axE1.set_xticks(pos)
axE1.set_xticklabels(labels)
#axE1.bar(labels, [dataPHYS["SpikeDurAvg"], eMSE_ISO, spikeMSE_NG, spikeMSE_WAVE], color=[:black, :red, :blue, :green])
ylabel("Duration (s)")

jitter_dataPHYS_burst = 0.5 .* rand(length(phys_burst_durs))
jitter_dataISO_burst = 0.5 .* rand(length(dataISO["BurstDurs"])) .+ 1
jitter_dataEC_burst = 0.5 .* rand(length(dataEC["BurstDurs"])) .+ 2
jitter_dataWAVE_burst = 0.5 .* rand(length(dataWAVE["BurstDurs"])) .+ 3

axE2 = fig4.add_subplot(gs[5, 2])
axE2.scatter(jitter_dataPHYS_burst, phys_burst_durs ./ 1000, s=5.0, c=:black)
axE2.scatter(jitter_dataISO_burst, dataISO["BurstDurs"] ./ 1000, s=5.0, c=:red)
axE2.scatter(jitter_dataEC_burst, dataEC["BurstDurs"] ./ 1000, s=5.0, c=:blue)
axE2.scatter(jitter_dataWAVE_burst, dataWAVE["BurstDurs"] ./ 1000, s=5.0, c=:green)
axE2.set_xticks(pos)
axE2.set_xticklabels(labels)
#axE2.bar(labels, [burstMSE_ISO, burstMSE_NG, burstMSE_WAVE], color=[:red, :blue, :green])

jitter_dataPHYS_ibi = 0.5 .* rand(length(phys_ibis))
jitter_dataISO_ibi = 0.5 .* rand(length(dataISO["IBIs"])) .+ 1
jitter_dataEC_ibi = 0.5 .* rand(length(dataEC["IBIs"])) .+ 2
jitter_dataWAVE_ibi = 0.5 .* rand(length(dataWAVE["IBIs"])) .+ 3

axE3 = fig4.add_subplot(gs[5, 3])
axE3.scatter(jitter_dataPHYS_ibi, phys_ibis ./ 1000, s=5.0, c=:black)
axE3.scatter(jitter_dataISO_ibi, dataISO["IBIs"] ./ 1000, s=5.0, c=:red)
axE3.scatter(jitter_dataEC_ibi, dataEC["IBIs"] ./ 1000, s=5.0, c=:blue)
axE3.scatter(jitter_dataWAVE_ibi, dataWAVE["IBIs"] ./ 1000, s=5.0, c=:green)
axE3.set_xticks(pos)
axE3.set_xticklabels(labels)
#axE3.bar(labels, [IBIMSE_ISO, IBIMSE_NG, IBIMSE_WAVE], color=[:red, :blue, :green])

axE3.annotate("A", (0.01, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("B", (0.01, 0.75), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("C", (0.01, 0.60), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("D", (0.01, 0.40), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("E", (0.01, 0.20), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

#%%
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
#loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2021 A Computational Model - Sci. Rep\Figures"
print("[$(now())]: Saving the figure 4...")
fig4.savefig("$(loc)/Figure6_ModelComparison.png")
plt.close("all")
println(" Completed")

#%% Print the results of all the traces
println("Data Values")
println("-----------------------------------------------")
println("""

ISO_n = $(length(iso_baselines))
ISO Baseline Avg: $(iso_baseline_avg|>x->round(x, digits = 2))±$(iso_baseline_SEM|>x->round(x, digits = 2))
ISO Minimum Avg: $(iso_min_avg|>x->round(x, digits = 2))±$(iso_min_SEM|>x->round(x, digits = 2))
ISO Maximum Avg: $(iso_max_avg|>x->round(x, digits = 2))±$(iso_max_SEM|>x->round(x, digits = 2))
ISO Spike Duration Avg: $(dataISO["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataISO["SpikeDurSEM"]|>x->round(x, digits = 2))
ISO Burst Duration Avg: $(dataISO["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataISO["BurstDurSEM"]|>x->round(x, digits = 2))
ISO IBI_Avg: $(dataISO["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataISO["IBISEM"]./1000|>x->round(x, digits = 2))

NO GABA n = $(length(ng_baselines))
NO GABA Baseline Avg: $(ng_baseline_avg|>x->round(x, digits = 2))±$(ng_baseline_SEM|>x->round(x, digits = 2))
NO GABA Minimum Avg: $(ng_min_avg|>x->round(x, digits = 2))±$(ng_min_SEM|>x->round(x, digits = 2))
NO GABA Maximum Avg: $(ng_max_avg|>x->round(x, digits = 2))±$(ng_max_SEM|>x->round(x, digits = 2))
NO GABA Spike Duration Avg: $(dataNG["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataNG["SpikeDurSEM"]|>x->round(x, digits = 2))
NO GABA Burst Duration Avg: $(dataNG["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataNG["BurstDurSEM"]|>x->round(x, digits = 2))
NO GABA IBI_Avg: $(dataNG["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataNG["IBISEM"]./1000|>x->round(x, digits = 2))

ECl -55 n = $(length(ec_baselines))
ECl -55 Baseline Avg: $(ec_baseline_avg|>x->round(x, digits = 2))±$(ec_baseline_SEM|>x->round(x, digits = 2))
ECl -55 Minimum Avg: $(ec_min_avg|>x->round(x, digits = 2))±$(ec_min_SEM|>x->round(x, digits = 2))
ECl -55 Maximum Avg: $(ec_max_avg|>x->round(x, digits = 2))±$(ec_max_SEM|>x->round(x, digits = 2))
ECl -55 Spike Duration Avg: $(dataEC["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataEC["SpikeDurSEM"]|>x->round(x, digits = 2))
ECl -55 Burst Duration Avg: $(dataEC["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataEC["BurstDurSEM"]|>x->round(x, digits = 2))
ECl -55 IBI_Avg: $(dataEC["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataEC["IBISEM"]./1000|>x->round(x, digits = 2))

WAVE n = $(length(wave_baselines))
WAVE Baseline Avg: $(wave_baseline_avg|>x->round(x, digits = 2))±$(wave_baseline_SEM|>x->round(x, digits = 2))
WAVE Minimum Avg: $(wave_min_avg|>x->round(x, digits = 2))±$(wave_min_SEM|>x->round(x, digits = 2))
WAVE Maximum Avg: $(wave_max_avg|>x->round(x, digits = 2))±$(wave_max_SEM|>x->round(x, digits = 2))
WAVE Spike Duration Avg: $(dataWAVE["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataWAVE["SpikeDurSEM"]|>x->round(x, digits = 2))
WAVE Burst Duration Avg: $(dataWAVE["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataWAVE["BurstDurSEM"]|>x->round(x, digits = 2))
WAVE IBI_Avg: $(dataWAVE["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataWAVE["IBISEM"]./1000|>x->round(x, digits = 2))

PHYS n = $(length(phys_baseline))
PHYS Baseline Avg: $(phys_baseline_avg|>x->round(x, digits = 2))±$(phys_baseline_SEM|>x->round(x, digits = 2))
PHYS Minimum Avg: $(phys_min_avg|>x->round(x, digits = 2))±$(phys_min_SEM|>x->round(x, digits = 2))
PHYS Maximum Avg: $(phys_max_avg|>x->round(x, digits = 2))±$(phys_max_SEM|>x->round(x, digits = 2))
PHYS Spike Duration Avg: $(phys_spike_avg|>x->round(x, digits = 2))±$(phys_spike_SEM |>x->round(x, digits = 2))
PHYS Burst Duration Avg: $(phys_burst_avg|>x->round(x, digits = 2))±$(phys_burst_SEM|>x->round(x, digits = 2))
PHYS IBI_Avg: $(phys_ibis_avg./1000|>x->round(x, digits = 2))±$(phys_ibis_SEM./1000|>x->round(x, digits = 2))
""")