using Revise
using RetinalChaos
#load the figure data from figure 5
include("figure_setup.jl")
include("making_figure5.jl")

#%% Set up the animation
maximum(ACH_map)

animate_dt = 50
animate_stops = 1:animate_dt:size(ACH_map, 3)
n_frames = length(animate_stops)
# Initialize the image plot
#All videos should be 4:3, 16:9, or 21:9
im1_levels = 0.0:0.05:2.0
im2_levels = 0.0:0.05:2.0
im3_levels = -13.0:0.5:13.0
figsize = (21, 9) #These are common video formats
fig, ax = subplots(1, 3, figsize=figsize)
plt.subplots_adjust(right=0.91, left=0.1, top=0.95, bottom=0.08, wspace=0.20)

ax[1].annotate("A", (0.07, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=25.0, fontweight="bold")
ax[2].annotate("B", (0.37, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=25.0, fontweight="bold")
ax[3].annotate("C", (0.65, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=25.0, fontweight="bold")

im1 = ax[1].contourf(ACH_map[:, :, 1], levels=im1_levels, vmin=0.0, vmax=0.75, cmap="viridis", extend="max")
ax[1].set_aspect("equal")
ax[1].set_title("Acetylcholine", fontsize=20.0)
cbare = fig.colorbar(im1, ax=ax[1], ticks=[0.0, 0.50, 1.0], aspect=5, location="bottom", pad=0.15)
cbare.ax.set_xlabel("ACh Release (mM)")
ax[1].spines["bottom"].set_visible(false)
ax[1].spines["left"].set_visible(false)
ax[1].set_xlabel("Cell # x", fontsize=15.0)
ax[1].set_ylabel("Cell # y", fontsize=15.0)

im2 = ax[2].contourf(GABA_map[:, :, 1], levels=im2_levels, vmin=0.0, vmax=0.75, cmap="jet", extend="max")
ax[2].set_aspect("equal")
cbari = fig.colorbar(im2, ax=ax[2], ticks=[0.0, 0.50, 1.0], aspect=5, location="bottom", pad=0.15)
cbari.ax.set_xlabel("GABA Release (mM)")
ax[2].set_title("GABA", fontsize=20.0)
ax[2].spines["bottom"].set_visible(false)
ax[2].spines["left"].set_visible(false)
ax[2].set_xlabel("Cell # x", fontsize=15.0)

im3 = ax[3].contourf(I_TOTAL[:, :, 1], levels=im3_levels, vmin=-5.0, vmax=5.0, cmap="RdYlGn", extend="both")
ax[3].set_aspect("equal")
cbarI = fig.colorbar(im3, ax=ax[3], ticks=[-4.0, 0.0, 4.0], aspect=5, location="bottom", pad=0.15)
cbarI.ax.set_xlabel("Induced Current (pA) ")
ax[3].set_title("Induced Current", fontsize=20.0)
ax[3].spines["bottom"].set_visible(false)
ax[3].spines["left"].set_visible(false)
ax[3].set_xlabel("Cell # x", fontsize=15.0)
ann1 = ax[1].annotate("t = 0.0 s", (0.05, 0.25), xycoords="figure fraction", annotation_clip=false, fontsize=25.0, fontweight="bold")

function animateinit() #tells our animator what artists will need re-drawing every time
     return im1, im2, im3, ann1
end
# Function to update the plot for each frame of the animation
function update(idx)
     print("Animating frame $idx of $n_frames ")
     print(size(ax[1].collections))
     print(", ")
     print(size(ax[2].collections))
     print(", ")
     println(size(ax[3].collections))

     for tp in ax[1].collections
          tp.remove()
     end
     for tp in ax[2].collections
          tp.remove()
     end
     for tp in ax[3].collections
          tp.remove()
     end

     tstop = animate_stops[idx+1]
     im1 = ax[1].contourf(ACH_map[:, :, tstop], levels=im1_levels, cmap="viridis")
     im2 = ax[2].contourf(GABA_map[:, :, tstop], levels=im2_levels, cmap="jet")
     im3 = ax[3].contourf(I_TOTAL[:, :, tstop], levels=im3_levels, vmin=-6.0, vmax=6.0, cmap="RdYlGn")
     ann1.set_text("t = $(round(tstop/1000, digits = 2)) s")
     return im1, im2, im3, ann1
end

# Set up the animation using PyPlot's animation module
animation = anim.FuncAnimation(fig, update, init_func=animateinit, frames=n_frames, repeat=false)
#plt.show()
# Save the animation as a gif using imagemagick
animation.save("$(loc)\\Supplementary Movie 1 Neurotransmitter Conductance.mp4", fps=1000.0 / animate_dt, dpi=100, writer="ffmpeg")
#%%
plt.close("all")

#%% Animate the Wave models
using RetinalChaos
include("figure_setup.jl")

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling"
save_path = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 2\Figures"

isolated_path = "$(data_root)/isolated_model"
dataISO = load("$(isolated_path)/data.jld2")
aniISO = animate_solution(dataISO, save_path; title = "All neurotransmission blocked", animation_filename="Supplementary Movie 2 All neurotransmission blocked")
plt.close("all")

noGABA_path = "$(data_root)/no_GABA_model"
dataNG = load("$(noGABA_path)/data.jld2")
animate_solution(dataNG, save_path; title = "No GABA transmission", animation_filename="Supplementary Movie 3 No GABA transmission")
plt.close("all")

wave_path = "$(data_root)\\wave_model"
dataHYP = load("$(wave_path)/data.jld2")
animate_solution(dataHYP, save_path; title="ECl = -65mV", animation_filename="Supplementary Movie 4 Hyperpolarizing GABA channel simulation")
plt.close("all")

ec_path = "$(data_root)\\ECl55_model"
dataDEP = load("$(ec_path)/data.jld2")
animate_solution(dataDEP, save_path; title="ECl = -55mV", animation_filename="Supplementary Movie 5 Depolarizing GABA channel simulation")
plt.close("all")