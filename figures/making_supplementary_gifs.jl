using Revise

#load the figure data from figure 5
include("making_figure5.jl")

#%% Set up the animation
maximum(GABA_map)
animate_dt = 50
animate_stops = 1:animate_dt:size(ACH_map, 3)
n_frames = length(animate_stops)

# Initialize the image plot
figsize = (7.5, 4)
fig, ax = subplots(1, 3, figsize=figsize)
plt.subplots_adjust(right=0.91, left=0.1, top=0.95, bottom=0.08, wspace=0.20)

ax[1].annotate("A", (0.07, 0.85), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
ax[2].annotate("B", (0.37, 0.85), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
ax[3].annotate("C", (0.65, 0.85), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
im1 = ax[1].contourf(ACH_map[:, :, 1], levels=0.0:0.05:1.0, cmap="viridis", extend="max")
ax[1].set_aspect("equal")
ax[1].set_title("Acetylcholine")
cbare = fig.colorbar(im1, ax=ax[1], ticks=[0.0, 0.25, 0.75], aspect=5, location="bottom", pad=0.15)
cbare.ax.set_xlabel("ACh Release (mM)")
ax[1].spines["bottom"].set_visible(false)
ax[1].spines["left"].set_visible(false)
ax[1].set_xlabel("Cell # x")
ax[1].set_ylabel("Cell # y")

im2 = ax[2].contourf(GABA_map[:, :, 1], levels=0.0:0.05:1.0, cmap="jet", extend="max")
ax[2].set_aspect("equal")
cbari = fig.colorbar(im2, ax=ax[2], ticks=[0.0, 0.25, 0.75], aspect=5, location="bottom", pad=0.15)
cbari.ax.set_xlabel("GABA Release (mM)")
ax[2].set_title("GABA")
ax[2].spines["bottom"].set_visible(false)
ax[2].spines["left"].set_visible(false)
ax[2].set_xlabel("Cell # x")

im3 = ax[3].contourf(I_TOTAL[:, :, 1], levels=-6.0:0.5:6.0, cmap="RdYlGn", extend="both")
ax[3].set_aspect("equal")
cbarI = fig.colorbar(im3, ax=ax[3], ticks=[-4.0, 0.0, 4.0], aspect=5, location="bottom", pad=0.15)
cbarI.ax.set_xlabel("Induced Current (pA) ")
ax[3].set_title("Induced Current")
ax[3].spines["bottom"].set_visible(false)
ax[3].spines["left"].set_visible(false)
ax[3].set_xlabel("Cell # x")
ann1 = ax[1].annotate("", (0.02, 0.24), xycoords="figure fraction", annotation_clip=false, fontweight="bold")

# Function to update the plot for each frame of the animation
function update(idx)
     println("Animating frame $idx of $n_frames")
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
     e_frame = ACH_map[:, :, tstop]
     i_frame = GABA_map[:, :, tstop]
     I_frame = I_TOTAL[:, :, tstop]
     im1 = ax[1].contourf(e_frame, levels=0.0:0.01:0.5, cmap="viridis")
     im2 = ax[2].contourf(i_frame, levels=0.0:0.01:0.5, cmap="jet")
     im3 = ax[3].contourf(I_frame, levels=-10.0:0.5:10.0, cmap="RdYlGn")
     ann1.set_text("t = $(round(tstop/1000, digits = 2)) s")
     return im1, im2, im3, ann1
end
# Set up the animation using PyPlot's animation module
animation = anim.FuncAnimation(fig, update, frames=n_frames, repeat=false)
# Save the animation as a gif using imagemagick
animation.save("$(loc)/Supplementary GIF 1 Neurotransmitter Conductance.gif", fps=1000.0 / animate_dt, dpi=100, writer="imagemagick")

#%% Animate the no 
figsize = (7.5, 4)
fig, ax = subplots(1, 3, figsize=figsize)
