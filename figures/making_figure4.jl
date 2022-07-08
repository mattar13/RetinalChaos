#=
=#
using Revise
using RetinalChaos
using PhysAnalysis, ABFReader
include("figure_setup.jl")
include("opening_data.jl")
#%%We have to open 3 different datasets

#1) open the physiological data loaded from the single cell recordings I made
loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
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

gs = fig4.add_gridspec(5, 2,
     width_ratios=(0.50, 0.50),
     height_ratios=(0.15, 0.15, 0.15, 0.275, 0.275),
     right=0.95, left=0.07,
     top=0.95, bottom=0.08,
     wspace=0.10, hspace=0.4)

axA = fig4.add_subplot(py"""$(gs)[0, :]""")
axA.plot(t_phys, data.data_array[1, idxs, 1], c=:green)
axA.xaxis.set_visible(false) #Turn off the bottom axis
axA.spines["bottom"].set_visible(false)

axB = fig4.add_subplot(py"""$(gs)[1, :]""")
axB.plot(t_iso, vt_iso, c=:blue)
axB.xaxis.set_visible(false) #Turn off the bottom axis
axB.spines["bottom"].set_visible(false)

axC = fig4.add_subplot(py"""$(gs)[2, :]""")
axC.plot(t_reg, vt_reg, c=:red)

#Figure D spike duration
gsDR = gs[4, 1].subgridspec(ncols=1, nrows=3)
axDR1 = fig4.add_subplot(gsDR[1, 1])
axDR1.plot(sdur_edges, sdur_weights, c=:green)
axDR1.xaxis.set_visible(false) #Turn off the bottom axis
axDR1.spines["bottom"].set_visible(false)

axDR2 = fig4.add_subplot(gsDR[2, 1])
axDR2.plot(iso_sdur_edges, iso_sdur_weights, c=:blue)
axDR2.xaxis.set_visible(false) #Turn off the bottom axis
axDR2.spines["bottom"].set_visible(false)

axDR3 = fig4.add_subplot(gsDR[3, 1])
axDR3.plot(wave_sdur_edges, wave_sdur_weights, c=:red)

gsDL = gs[4, 2].subgridspec(ncols=1, nrows=3)

axDL1 = fig4.add_subplot(gsDL[1, 1])
axDL1.plot(isi_edges, isi_weights, c=:green)
axDL1.xaxis.set_visible(false) #Turn off the bottom axis
axDL1.spines["bottom"].set_visible(false)

axDL2 = fig4.add_subplot(gsDL[2, 1])
axDL2.plot(iso_isi_edges, iso_isi_weights, c=:blue)
axDL2.xaxis.set_visible(false) #Turn off the bottom axis
axDL2.spines["bottom"].set_visible(false)

axDL3 = fig4.add_subplot(gsDL[3, 1])
axDL3.plot(wave_isi_edges, wave_isi_weights, c=:red)

gsER = gs[5, 1].subgridspec(ncols=1, nrows=3)

axER1 = fig4.add_subplot(gsER[1, 1])
axER1.plot(bdur_edges ./ 1000, bdur_weights, c=:green)

axER2 = fig4.add_subplot(gsER[2, 1])
axER2.plot(iso_bdur_edges ./ 1000, iso_bdur_weights, c=:blue)

axER3 = fig4.add_subplot(gsER[3, 1])
axER3.plot(wave_bdur_edges ./ 1000, wave_bdur_weights, c=:red)

gsEL = gs[5, 2].subgridspec(ncols=1, nrows=3)

axEL1 = fig4.add_subplot(gsEL[1, 1])
axEL1.plot(ibi_edges, ibi_weights, c=:green)

axEL2 = fig4.add_subplot(gsEL[2, 1])
axEL2.plot(iso_ibi_edges, iso_ibi_weights, c=:blue)

axEL3 = fig4.add_subplot(gsEL[3, 1])
axEL3.plot(wave_ibi_edges, wave_ibi_weights, c=:red)