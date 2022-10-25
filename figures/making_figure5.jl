using Revise
include("figure_setup.jl")
using RetinalChaos
using ePhys
include("opening_data.jl")

#%%
iso_t = isolated_data["Time"]
iso_raster = isolated_data["DataArray"]
iso_spikes = iso_raster .> isolated_data["Thresholds"]
iso_vt = reshape(isolated_data["DataArray"], 64, 64, size(isolated_data["DataArray"], 2))
averaged_iso_vt = sum(iso_vt, dims=3) / size(iso_vt, 3)
counted_iso = iso_spikes' .* collect(1:size(iso_spikes, 2))
counted_iso = reshape(counted_iso', 64, 64, size(iso_spikes, 2))
highest_iso = maximum(counted_iso, dims=3) ./ 1000


noGABA_t = noGABA_data["Time"]
noGABA_raster = noGABA_data["DataArray"]
noGABA_spikes = noGABA_raster .> noGABA_data["Thresholds"]
noGABA_vt = reshape(noGABA_data["DataArray"], 64, 64, size(noGABA_data["DataArray"], 2))
averaged_noGABA_vt = sum(noGABA_vt, dims=3) / size(noGABA_vt, 3)
counted_noGABA = noGABA_spikes .* collect(1:size(noGABA_spikes, 2))'
counted_noGABA = reshape(counted_noGABA[:, 1:20000], 64, 64, 20000)
highest_noGABA = maximum(counted_noGABA, dims=3) ./ 1000

wave_t = wave_data["Time"]
wave_raster = wave_data["DataArray"]
wave_spikes = wave_raster .> wave_data["Thresholds"]
wave_vt = reshape(wave_data["DataArray"], 64, 64, size(wave_data["DataArray"], 2))
averaged_wave_vt = sum(wave_vt, dims=3) / size(wave_vt, 3)
counted_wave = wave_spikes .* collect(1:size(wave_spikes, 2))'
counted_wave = reshape(counted_wave[:, 20001:40000], 64, 64, 20000)
highest_wave = maximum(counted_wave, dims=3) ./ 1000


function make_heatmap_plot(axes, data, idx; levels=-95.0:10.0:5, xlabels = nothing)
     axes.set_aspect("equal", "box")
     axes.yaxis.set_visible(false) #Turn off the bottom axis
     axes.spines["right"].set_visible(false)
     if isnothing(xlabel)
          axes.xaxis.set_visible(false) #Turn off the bottom axis
          axes.spines["bottom"].set_visible(false)
     else
          #axes.xaxis.set_visible(false) #Turn off the bottom axis
          #axes.spines["bottom"].set_visible(false)
          axes.set_xticks([])
          xlabel(xlabels)
     end
     ctr = axes.contourf(data[:, :, idx], levels=levels, extend="both")
     return ctr
end

#%% We are going to make figure 5 here
print("[$(now())]: Plotting... ")
width_inches = 15.0
height_inches = 18.0
fig5 = plt.figure("Wave Data", figsize=(width_inches, height_inches))


gs = fig5.add_gridspec(6, 6,
     #width_ratios = (0.18, 0.18, 0.18, 0.18, 0.18, 0.10), 
     right=0.95, left=0.10,
     top=0.95, bottom=0.08,
     wspace=0.38, hspace=0.6)


axA = fig5.add_subplot(py"""$(gs)[0, :]""")
#pick 10 random indexes to plot
trace_idxs = rand(1:size(iso_vt, 1), (10))
for idx in trace_idxs
     axA.plot(iso_t ./ 1000, iso_raster[idx, :], lw=3.0)
end
ylabel("Vt (mv)")
xlabel("Time (s)")

axB1 = fig5.add_subplot(gs[2, 1])
make_heatmap_plot(axB1, iso_vt, 4000, xlabels="t=4s")

axB2 = fig5.add_subplot(gs[2, 2])
make_heatmap_plot(axB2, iso_vt, 6000, xlabels = "t=62")

axB3 = fig5.add_subplot(gs[2, 3])
make_heatmap_plot(axB3, iso_vt, 8000, xlabels = "t=8s")

axB4 = fig5.add_subplot(gs[2, 4])
make_heatmap_plot(axB4, iso_vt, 10000, xlabels = "t=10s")

axB5 = fig5.add_subplot(gs[2, 5])
ctrB5 = make_heatmap_plot(axB5, iso_vt, 12000, xlabels="t=12s")
cbarB5 = fig5.colorbar(ctrB5, ticks=collect(-95.0:25.0:5))
cbarB5.ax.set_ylabel("Voltage (mV)")

axB6 = fig5.add_subplot(gs[2, 6])
ctrB6 = make_heatmap_plot(axB6, highest_iso, 1; levels=1:10:120)
cbarB6 = fig5.colorbar(ctrB6, ticks=0:30:120)
cbarB6.ax.set_ylabel("Time (s)")

axC = fig5.add_subplot(py"""$(gs)[2, :]""")
for idx in trace_idxs
     axC.plot(noGABA_t ./ 1000, noGABA_raster[idx, :], lw=3.0)
end
ylabel("Vt (mv)")
xlabel("Time (s)")

axD1 = fig5.add_subplot(gs[4, 1])
make_heatmap_plot(axD1, noGABA_vt, 10000, xlabels = "t=10s")

axD2 = fig5.add_subplot(gs[4, 2])
make_heatmap_plot(axD2, noGABA_vt, 12000, xlabels = "t=12s")

axD3 = fig5.add_subplot(gs[4, 3])
make_heatmap_plot(axD3, noGABA_vt, 14000, xlabels = "t=14s")

axD4 = fig5.add_subplot(gs[4, 4])
make_heatmap_plot(axD4, noGABA_vt, 16000, xlabels = "t=16s")

axD5 = fig5.add_subplot(gs[4, 5])
ctrD5 = make_heatmap_plot(axD5, noGABA_vt, 18000, xlabels="t=16s")
cbarD5 = fig5.colorbar(ctrD5, ticks=collect(-95.0:45.0:5))
cbarD5.ax.set_ylabel("Voltage (mV)")

axD6 = fig5.add_subplot(gs[4, 6])
ctrD6 = make_heatmap_plot(axD6, highest_noGABA, 1; levels=1:10:120)
cbarD6 = fig5.colorbar(ctrD6, ticks=0:30:120)
cbarD6.ax.set_ylabel("Time (s)")

axE = fig5.add_subplot(py"""$(gs)[4, :]""")
for idx in trace_idxs
     axE.plot(wave_t ./ 1000, wave_raster[idx, :], lw=3.0)
end
ylabel("Vt (mv)")
xlabel("Time (s)")

axF1 = fig5.add_subplot(gs[6, 1])
make_heatmap_plot(axF1, wave_vt, 20000, xlabels = "t=20s")

axF2 = fig5.add_subplot(gs[6, 2])
make_heatmap_plot(axF2, wave_vt, 22000, xlabels = "t=22s")

axF3 = fig5.add_subplot(gs[6, 3])
make_heatmap_plot(axF3, wave_vt, 24000, xlabels = "t=24s")

axF4 = fig5.add_subplot(gs[6, 4])
make_heatmap_plot(axF4, wave_vt, 26000, xlabels = "t=26s")

axF5 = fig5.add_subplot(gs[6, 5])
ctrF5 = make_heatmap_plot(axF5, wave_vt, 28000, xlabels="t=28s")
cbarF5 = fig5.colorbar(ctrF5, ticks=collect(-95.0:45.0:5))
cbarF5.ax.set_ylabel("Voltage (mV)")

axF6 = fig5.add_subplot(gs[6, 6])
ctrF6 = make_heatmap_plot(axF6, highest_wave, 1; levels=1:10:120)
cbarF6 = fig5.colorbar(ctrF6, ticks=0:30:120)
cbarF6.ax.set_ylabel("Time (s)")

axF6.annotate("A", (0.01, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axF6.annotate("B", (0.01, 0.80), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axF6.annotate("C", (0.01, 0.65), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axF6.annotate("D", (0.01, 0.50), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axF6.annotate("E", (0.01, 0.35), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axF6.annotate("F", (0.01, 0.20), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")

#%% Shit
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2021 A Computational Model - Sci. Rep\Figures"
print("[$(now())]: Saving the figure 5...")
fig5.savefig("$(loc)/figure5_WaveModel.png")
plt.close("All")
println(" Completed")