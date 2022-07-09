#=
=#
using Revise
using RetinalChaos
using PhysAnalysis, ABFReader
include("figure_setup.jl")
include("opening_data.jl")
#%%We have to open 3 different datasets

#1) open the physiological data loaded from the single cell recordings I made
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
data = readABF(target_file, channels=["Vm_prime4"], stimulus_name=nothing, time_unit=:ms)
data - 25.0
timestamps, _ = timeseries_analysis(data.t, data.data_array[:, :, 1])
bursts = timestamps["Bursts"][1]

t_phys = bursts[2, 1]-1000:1.0:bursts[2, 1]+119e3
idxs = round.(Int64, t_phys ./ data.dt)
t_phys = t_phys .- t_phys[1]

#2) Open the isolated
iso_arr = isolated_data["DataArray"]
t_iso = 1:size(iso_arr, 2)
vt_iso = iso_arr[rand(1:size(iso_arr, 1)), :]

#3) Open the wave model
reg_arr = wave_data["DataArray"]
t_reg = 1:size(reg_arr, 2)
vt_reg = reg_arr[rand(1:size(reg_arr, 1)), :]
t_reg = t_reg ./ 1000.0

#%%
print("[$(now())]: Plotting... ")
width_inches = 15.0
height_inches = 10.0
fig4 = plt.figure("Physiology Data", figsize=(width_inches, height_inches))

col1_ylabel = -0.05
col2_ylabel = -0.11

gs = fig4.add_gridspec(5, 2,
     width_ratios=(0.50, 0.50),
     height_ratios=(0.15, 0.15, 0.15, 0.275, 0.275),
     right=0.95, left=0.10,
     top=0.95, bottom=0.08,
     wspace=0.10, hspace=0.6)

axA = fig4.add_subplot(py"""$(gs)[0, :]""")
ylabel("Vm (mV)")
axA.plot(t_phys, data.data_array[1, idxs, 1], c=:green, lw=3.0)
axA.xaxis.set_visible(false) #Turn off the bottom axis
axA.spines["bottom"].set_visible(false)
axA.yaxis.set_label_coords(col1_ylabel, 0.5)

axB = fig4.add_subplot(py"""$(gs)[1, :]""")
ylabel("Vt (mV)")
axB.plot(t_iso, vt_iso, c=:blue, lw=3.0)
axB.xaxis.set_visible(false) #Turn off the bottom axis
axB.spines["bottom"].set_visible(false)
axB.yaxis.set_label_coords(col1_ylabel, 0.5)

axC = fig4.add_subplot(py"""$(gs)[2, :]""")
ylabel("Vt (mV)")
axC.plot(t_reg, vt_reg, c=:red, lw=3.0)
axC.yaxis.set_label_coords(col1_ylabel, 0.5)

#Figure D spike duration
gsDR = gs[4, 1].subgridspec(ncols=1, nrows=3)

axDR1 = fig4.add_subplot(gsDR[1, 1])
axDR1.plot(sdur_edges, sdur_weights, c=:green, lw=3.0)
axDR1.xaxis.set_visible(false) #Turn off the bottom axis
axDR1.spines["bottom"].set_visible(false)

axDR2 = fig4.add_subplot(gsDR[2, 1])
ylabel("Probability")
axDR2.plot(iso_sdur_edges, iso_sdur_weights, c=:blue, lw=3.0)
axDR2.xaxis.set_visible(false) #Turn off the bottom axis
axDR2.spines["bottom"].set_visible(false)
axDR2.yaxis.set_label_coords(col2_ylabel, 0.5)

axDR3 = fig4.add_subplot(gsDR[3, 1])
axDR3.plot(wave_sdur_edges, wave_sdur_weights, c=:red, lw=3.0)

gsDL = gs[4, 2].subgridspec(ncols=1, nrows=3)

axDL1 = fig4.add_subplot(gsDL[1, 1])
axDL1.plot(isi_edges, isi_weights, c=:green, lw=3.0)
axDL1.xaxis.set_visible(false) #Turn off the bottom axis
axDL1.spines["bottom"].set_visible(false)

axDL2 = fig4.add_subplot(gsDL[2, 1])
axDL2.plot(iso_isi_edges, iso_isi_weights, c=:blue, lw=3.0)
axDL2.xaxis.set_visible(false) #Turn off the bottom axis
axDL2.spines["bottom"].set_visible(false)

axDL3 = fig4.add_subplot(gsDL[3, 1])
axDL3.plot(wave_isi_edges, wave_isi_weights, c=:red, lw=3.0)

gsER = gs[5, 1].subgridspec(ncols=1, nrows=3)

axER1 = fig4.add_subplot(gsER[1, 1])
axER1.plot(bdur_edges ./ 1000, bdur_weights, c=:green, lw=3.0)
axER1.xaxis.set_visible(false) #Turn off the bottom axis
axER1.spines["bottom"].set_visible(false)

axER2 = fig4.add_subplot(gsER[2, 1])
ylabel("Probability")
axER2.plot(iso_bdur_edges ./ 1000, iso_bdur_weights, c=:blue, lw=3.0)
axER2.xaxis.set_visible(false) #Turn off the bottom axis
axER2.spines["bottom"].set_visible(false)
axER2.yaxis.set_label_coords(col2_ylabel, 0.5)

axER3 = fig4.add_subplot(gsER[3, 1])
axER3.plot(wave_bdur_edges ./ 1000, wave_bdur_weights, c=:red, lw=3.0)

gsEL = gs[5, 2].subgridspec(ncols=1, nrows=3)

axEL1 = fig4.add_subplot(gsEL[1, 1])
axEL1.plot(ibi_edges / 1000, ibi_weights, c=:green, lw=3.0)
axEL1.xaxis.set_visible(false) #Turn off the bottom axis
axEL1.spines["bottom"].set_visible(false)

axEL2 = fig4.add_subplot(gsEL[2, 1])
axEL2.plot(iso_ibi_edges / 1000, iso_ibi_weights, c=:blue, lw=3.0)
axEL2.xaxis.set_visible(false) #Turn off the bottom axis
axEL2.spines["bottom"].set_visible(false)

axEL3 = fig4.add_subplot(gsEL[3, 1])
axEL3.plot(wave_ibi_edges / 1000, wave_ibi_weights, c=:red, lw=3.0)

axEL3.annotate("A", (0.01, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axEL3.annotate("B", (0.01, 0.80), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axEL3.annotate("C", (0.01, 0.70), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axEL3.annotate("D", (0.01, 0.50), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")
axEL3.annotate("E", (0.01, 0.25), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")


#%%
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2021 A Computational Model - Sci. Rep\Figures"
print("[$(now())]: Saving the figure 4...")
fig4.savefig("$(loc)/figure4_ModelComparison.png")
plt.close("All")
println(" Completed")