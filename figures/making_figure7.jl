using Revise
include("figure_setup.jl")
using RetinalChaos
using ePhys
#include("opening_data.jl")

#4x4 grid, with raster plots of all spiking activity 
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"

wave_path = "$(data_root)/wave_model"
dataWAVE = load("$(wave_path)/data.jld2")
tsWAVE = load("$(wave_path)/timestamps.jld2")
tsWAVE = convert(Dict{String,Vector{Matrix{Float64}}}, tsWAVE)
WAVE_t = dataWAVE["Time"][1000:61000]
WAVE_raster = dataWAVE["DataArray"][1:200, 1000:61000]
WAVE_sDurs = dataWAVE["SpikeDurs"]
WAVE_bDurs = dataWAVE["BurstDurs"]
WAVE_IBI = dataWAVE["IBIs"]

FASTei_path = "$(data_root)/FastI_Ediffusion"
dataFASTei = load("$(FASTei_path)/data.jld2")
tsFASTei = load("$(FASTei_path)/timestamps.jld2")
tsFASTei = convert(Dict{String,Vector{Matrix{Float64}}}, tsFASTei)
FASTei_t = dataFASTei["Time"][1000:61000]
FASTei_raster = dataFASTei["DataArray"][1:200, 1000:61000]
FASTei_sDurs = dataFASTei["SpikeDurs"]
FASTei_bDurs = dataFASTei["BurstDurs"]
FASTei_IBI = dataFASTei["IBIs"]

SLOWei_path = wave_path = "$(data_root)/SlowI_Ediffusion"
dataSLOWei = load("$(SLOWei_path)/data.jld2")
tsSLOWei = load("$(SLOWei_path)/timestamps.jld2")
tsSLOWei = convert(Dict{String,Vector{Matrix{Float64}}}, tsSLOWei)
SLOWei_t = dataSLOWei["Time"][1000:61000]
SLOWei_raster = dataSLOWei["DataArray"][1:200, 1000:61000]
SLOWei_sDurs = dataSLOWei["SpikeDurs"]
SLOWei_bDurs = dataSLOWei["BurstDurs"]
SLOWei_IBI = dataSLOWei["IBIs"]

#%% We are going to make figure 5 here
levels = -95.0:5.0:5
print("[$(now())]: Plotting figure 5...")
width_inches = 7.5
height_inches = 5.0
fig5 = plt.figure("Wave Dynamics", figsize=(width_inches, height_inches))

gs = fig5.add_gridspec(2, 4,
     width_ratios=(0.325, 0.325, 0.325, 0.04),
     right=0.87, left=0.12,
     top=0.90, bottom=0.08,
     wspace=0.45, hspace=0.3
)

axA1 = fig5.add_subplot(gs[1, 1])
axA1.set_title("Fast (D = 0.01)")
ctrFASTei = axA1.contourf((FASTei_t .- FASTei_t[1]) ./ 1000, 1:200, FASTei_raster, levels=levels, extend="both", c=:curl)
axA1.set_ylabel("Cell Number")
axA1.xaxis.set_major_locator(MultipleLocator(20.0))
axA1.xaxis.set_minor_locator(MultipleLocator(10.0))
axA1.set_xlabel("Time (s)")

axA2 = fig5.add_subplot(gs[1, 2])
axA2.set_title("Medium (D = 0.005)")
ctrWAVE = axA2.contourf((WAVE_t .- WAVE_t[1]) ./ 1000, 1:200, WAVE_raster, levels=levels, extend="both", c=:curl)
axA2.yaxis.set_visible(false) #Turn off the bottom axis
axA2.spines["left"].set_visible(false)
axA2.xaxis.set_major_locator(MultipleLocator(20.0))
axA2.xaxis.set_minor_locator(MultipleLocator(10.0))
axA2.set_xlabel("Time (s)")

axA3 = fig5.add_subplot(gs[1, 3])
axA3.set_title("Slow (D = 0.001)")
ctrSLOWei = axA3.contourf((SLOWei_t .- SLOWei_t[1]) ./ 1000, 1:200, SLOWei_raster, levels=levels, extend="both", c=:curl)
axA3.yaxis.set_visible(false) #Turn off the bottom axis
axA3.spines["left"].set_visible(false)
axA3.xaxis.set_major_locator(MultipleLocator(20.0))
axA3.xaxis.set_minor_locator(MultipleLocator(10.0))
axA3.set_xlabel("Time (s)")

axA4 = fig5.add_subplot(gs[1, 4])
cbarSE = fig5.colorbar(ctrSLOWei, cax=axA4, ticks=[-90.0, 5.0])
cbarSE.ax.set_ylabel("Voltage (mV)")
axA1.annotate("A", (0.01, 0.90), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")

labels = ["FAST", "MED", "SLOW"]
pos = [0.25, 1.25, 2.25]
jitter_dataFASTei_spike = 0.5 .* rand(length(FASTei_sDurs))
jitter_dataWAVE_spike = 0.5 .* rand(length(WAVE_sDurs)) .+ 1
jitter_dataSLOWei_spike = 0.5 .* rand(length(SLOWei_sDurs)) .+ 2


axB1 = fig5.add_subplot(gs[2, 1])
axB1.set_ylabel("Spike Durations (ms)")
axB1.scatter(jitter_dataFASTei_spike, FASTei_sDurs, s=5.0, c=:green)
axB1.scatter(jitter_dataWAVE_spike, WAVE_sDurs, s=5.0, c=:blue)
axB1.scatter(jitter_dataSLOWei_spike, SLOWei_sDurs, s=5.0, c=:red)
axB1.set_xticks(pos)
axB1.set_xticklabels(labels)

jitter_dataFASTei_burst = 0.5 .* rand(length(FASTei_bDurs))
jitter_dataWAVE_burst = 0.5 .* rand(length(WAVE_bDurs)) .+ 1
jitter_dataSLOWei_burst = 0.5 .* rand(length(SLOWei_bDurs)) .+ 2

axB2 = fig5.add_subplot(gs[2, 2])
axB2.set_ylabel("Burst Durations (s)")
axB2.scatter(jitter_dataFASTei_burst, FASTei_bDurs ./ 1000, s=5.0, c=:green)
axB2.scatter(jitter_dataWAVE_burst, WAVE_bDurs ./ 1000, s=5.0, c=:blue)
axB2.scatter(jitter_dataSLOWei_burst, SLOWei_bDurs ./ 1000, s=5.0, c=:red)
axB2.set_xticks(pos)
axB2.set_xticklabels(labels)

jitter_dataFASTei_IBI = 0.5 .* rand(length(FASTei_IBI))
jitter_dataWAVE_IBI = 0.5 .* rand(length(WAVE_IBI)) .+ 1
jitter_dataSLOWei_IBI = 0.5 .* rand(length(SLOWei_IBI)) .+ 2

axB2 = fig5.add_subplot(gs[2, 3])
axB2.set_ylabel("IBIs (s)")
col1_ylabel = -0.15
axB2.yaxis.set_label_coords(col1_ylabel, 0.5)
axB2.scatter(jitter_dataFASTei_IBI, FASTei_IBI ./ 1000, s=5.0, c=:green)
axB2.scatter(jitter_dataWAVE_IBI, WAVE_IBI ./ 1000, s=5.0, c=:blue)
axB2.scatter(jitter_dataSLOWei_IBI, SLOWei_IBI ./ 1000, s=5.0, c=:red)
axB2.set_xticks(pos)
axB2.set_xticklabels(labels)

axB1.annotate("B", (0.01, 0.50), xycoords="figure fraction", annotation_clip=false, fontsize=30.0, fontweight="bold")

#%%
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
print("[$(now())]: Saving the figure 5...")
fig5.savefig("$(loc)/Figure7_DiffusionSpeeds.jpg")
plt.close("all")
println(" Completed")