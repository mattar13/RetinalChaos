using Revise
include("figure_setup.jl")
using Dates
using RetinalChaos
using StatsBase, Statistics
using ePhys
using DataFrames
import ePhys: dwt_filter
#include("opening_data.jl")

print("Extract all data points")
spike_t_PHYS, spike_vt_PHYS = extract_spike_trace(tsPHYS, dataPHYS, idx=5, dt=dt, normalize=false)
burst_t_PHYS, burst_vt_PHYS = extract_burst_trace(tsPHYS, dataPHYS, idx=2, dt=dt, normalize=false)
IBI_t_PHYS, IBI_vt_PHYS = extract_IBI_trace(tsPHYS, dataPHYS, idx=2, dt=dt, normalize=false)

spike_t_WAVE, spike_vt_WAVE = extract_spike_trace(tsWAVE, dataWAVE, normalize=false, cell_n=wave_xIdx, idx=1)
burst_t_WAVE, burst_vt_WAVE = extract_burst_trace(tsWAVE, dataWAVE, normalize=false, cell_n=wave_xIdx, idx=1)
IBI_t_WAVE, IBI_vt_WAVE = extract_IBI_trace(tsWAVE, dataWAVE, normalize=false, cell_n=wave_xIdx, idx=1)

spike_t_ISO, spike_vt_ISO = extract_spike_trace(tsISO, dataISO, normalize=false, cell_n=iso_xIdx, idx=1)
burst_t_ISO, burst_vt_ISO = extract_burst_trace(tsISO, dataISO, normalize=false, cell_n=iso_xIdx, idx=1)
IBI_t_ISO, IBI_vt_ISO = extract_IBI_trace(tsISO, dataISO, normalize=false, cell_n=iso_xIdx, idx=1)

spike_t_NG, spike_vt_NG = extract_spike_trace(tsNG, dataNG, normalize=false, cell_n=ng_xIdx, idx=1)
burst_t_NG, burst_vt_NG = extract_burst_trace(tsNG, dataNG, normalize=false, cell_n=ng_xIdx, idx=1)
IBI_t_NG, IBI_vt_NG = extract_IBI_trace(tsNG, dataNG, normalize=false, cell_n=ng_xIdx, idx=1)

spike_t_EC, spike_vt_EC = extract_spike_trace(tsEC, dataEC, normalize=false, cell_n=ec_xIdx, idx=1)
burst_t_EC, burst_vt_EC = extract_burst_trace(tsEC, dataEC, normalize=false, cell_n=ec_xIdx, idx=1)
IBI_t_EC, IBI_vt_EC = extract_IBI_trace(tsEC, dataEC, normalize=false, cell_n=ec_xIdx, idx=1)
println("Complete")

print("Extract wave rasters")
WAVE_t = dataWAVE["Time"][1000:61000]
WAVE_raster = dataWAVE["DataArray"][1:200, 1000:61000]

ISO_t = dataISO["Time"][1000:61000]
ISO_raster = dataISO["DataArray"][1:200, 1000:61000]

EC_t = dataEC["Time"][1000:61000]
EC_raster = dataEC["DataArray"][1:200, 1000:61000]
println("Complete")

#%% Open the wave stats from excel

#Need to open wave stats data
wavesWAVE = DataFrame(XLSX.readtable("$(wave_path)\\wave_data.xlsx", "Waves"))
wavesISO = DataFrame(XLSX.readtable("$(isolated_path)\\wave_data.xlsx", "Waves"))
wavesECl = DataFrame(XLSX.readtable("$(ECl55_path)\\wave_data.xlsx", "Waves"))

#%%
print("[$(now())]: Plotting... ")
width_inches = 7.5
height_inches = 7.5
fig4 = plt.figure("Physiology Data", figsize=(width_inches, height_inches))

col1_ylabel = -0.31
col2_ylabel = -0.11
levels = -95.0:10.0:5.0

gs = fig4.add_gridspec(5, 4,
     #width_ratios=(0.50, 0.50),
     #height_ratios=(0.16, 0.16, 0.15, 0.15, 0.15, 0.15),
     right=0.95, left=0.15,
     top=0.95, bottom=0.03,
     wspace=0.30, hspace=0.5)

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

axB4 = fig4.add_subplot(gs[2, 4])
axB4.contourf((ISO_t .- ISO_t[1]) ./ 1000, 1:200, ISO_raster, levels=levels, extend="both", color=:curl)
#axB4.set_ylabel("Cell Number")
axB4.xaxis.set_visible(false) #Turn off the bottom axis
axB4.spines["bottom"].set_visible(false)
axB4.xaxis.set_major_locator(MultipleLocator(20.0))
axB4.xaxis.set_minor_locator(MultipleLocator(10.0))
#axB4.set_xlabel("Time (s)")

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

axC4 = fig4.add_subplot(gs[3, 4])
axC4.contourf((EC_t .- EC_t[1]) ./ 1000, 1:200, EC_raster, levels=levels, extend="both", color=:curl)
#axB4.set_ylabel("Cell Number")
axC4.xaxis.set_visible(false) #Turn off the bottom axis
axC4.spines["bottom"].set_visible(false)
axC4.xaxis.set_major_locator(MultipleLocator(20.0))
axC4.xaxis.set_minor_locator(MultipleLocator(10.0))

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

axD4 = fig4.add_subplot(gs[4, 4])
ctrWAVE = axD4.contourf((WAVE_t .- WAVE_t[1]) ./ 1000, 1:200, WAVE_raster, levels=levels, extend="both", color=:curl)
#axB4.set_ylabel("Cell Number")
#axD4.xaxis.set_visible(false) #Turn off the bottom axis
#axD4.spines["bottom"].set_visible(false)
axD4.set_xlabel("Time (s)")
axD4.xaxis.set_major_locator(MultipleLocator(20.0))
axD4.xaxis.set_minor_locator(MultipleLocator(10.0))

axA4 = fig4.add_subplot(gs[1, 4])
axA4.set_title("Model Rasters")
cbarSE = fig4.colorbar(ctrWAVE, cax=axA4, ticks=[-90.0, 5.0])
cbarSE.ax.set_ylabel("Voltage (mV)")
#axA1.annotate("A", (0.01, 0.90), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")

#Make violin Graphs showing the MSE loss between each RandomNumbers
labels = ["P", "I", "D", "H"]
colors = [:black, :red, :blue, :green]
pos = [1.0, 2.0, 3.0, 4.0]

axE1 = fig4.add_subplot(gs[5, 1])
datasetSPIKE = [
     phys_spike_durs./1000, 
     dataISO["SpikeDurs"]./1000, 
     dataEC["SpikeDurs"]./1000, 
     dataWAVE["SpikeDurs"]./1000
]
vpSPIKE = axE1.violinplot(datasetSPIKE, widths = 0.9, showmeans = true, showextrema = true)
for (vp_idx, vp) in enumerate(vpSPIKE["bodies"])
     vp.set_facecolor(colors[vp_idx])
     vp.set_edgecolor(:black)
     vp.set_linewidth(1.0)
     vp.set_alpha(0.95)
end
axE1.set_xticks(pos)
axE1.set_xticklabels(labels)
axE1.set_ylabel("Duration (s)")

axE2 = fig4.add_subplot(gs[5, 2])
datasetBURST = [
     phys_burst_durs ./ 1000,
     dataISO["BurstDurs"] ./ 1000,
     dataEC["BurstDurs"] ./ 1000, 
     dataWAVE["BurstDurs"] ./ 1000
]
vpBURST = axE2.violinplot(datasetBURST, widths = 0.9, showmeans = true, showextrema = true)
for (vp_idx, vp) in enumerate(vpBURST["bodies"])
     vp.set_facecolor(colors[vp_idx])
     vp.set_edgecolor(:black)
     vp.set_linewidth(1.0)
     vp.set_alpha(0.95)
end
axE2.set_xticks(pos)
axE2.set_xticklabels(labels)

axE3 = fig4.add_subplot(gs[5, 3])
datasetIBI = [
     phys_ibis ./ 1000,
     dataISO["IBIs"] ./ 1000,
     dataEC["IBIs"] ./ 1000,
     dataWAVE["IBIs"] ./ 1000,
]
vpIBI = axE3.violinplot(datasetIBI, widths = 0.9, showmeans = true, showextrema = true)
for (vp_idx, vp) in enumerate(vpIBI["bodies"])
     vp.set_facecolor(colors[vp_idx])
     vp.set_edgecolor(:black)
     vp.set_linewidth(1.0)
     vp.set_alpha(0.95)
end
axE3.set_xticks(pos)
axE3.set_xticklabels(labels)

datasetWAVE = [
     [0.0], #Phys data is empty
     wavesWAVE.WaveTime./1000, 
     wavesISO.WaveTime./1000, 
     wavesECl.WaveTime./1000
]

axE4 = fig4.add_subplot(gs[5, 4])
vpWAVE = axE4.violinplot(datasetWAVE, widths = 0.9, showmeans = true, showextrema = true)
for (vp_idx, vp) in enumerate(vpWAVE["bodies"])
     vp.set_facecolor(colors[vp_idx])
     vp.set_edgecolor(:black)
     vp.set_linewidth(1.0)
     vp.set_alpha(0.95)
end
axE4.set_xticks(pos)
axE4.set_xticklabels(labels)
axE4.set_ylim(-0.1, 10.0)
#axE4.set_ylabel("Wave Duration(s)")
#axE3.bar(labels, [IBIMSE_ISO, IBIMSE_NG, IBIMSE_WAVE], color=[:red, :blue, :green])

axE3.annotate("A", (0.01, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("B", (0.01, 0.75), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("C", (0.01, 0.60), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("D", (0.01, 0.40), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axE3.annotate("E", (0.01, 0.20), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

#%%
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 2\Figures"
#loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2021 A Computational Model - Sci. Rep\Figures"
print("[$(now())]: Saving the figure 4...")
fig4.savefig("$(loc)/Figure6_ModelComparison.png")
plt.close("all")
println(" Completed")
