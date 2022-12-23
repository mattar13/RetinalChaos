#=
Figure 3 will include 4 rows. 
Row 1 will be a voltage trace with Acetylcholine and GABA activations
Row 2 will be GABA and ACh directional selectitivy
Row 3 will be current maps
Row 4 will be current induced by ACh and GABA
=#

using Revise

using RetinalChaos
import RetinalChaos: Φ, ħ, ∇
import RetinalChaos: T_ODE_NT_Release, get_timestamps, max_interval_algorithim
using PyPlot
include("figure_setup.jl")
#%% Simulating model
println("Simulating model for giure 5")

import RetinalChaos.ODEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters

print("[$(now())]: Setting up modelling data... ")
parameters[I_app] = 15.0

tmin = 0.0
dt = 1.0; #Set the time differential
tmax = 300e3
#Step 4: set up the problem
prob = ODEProblem(ODEModel, u0, (tmin, tmax), parameters)
#Step 5: Solve the problem
sol = solve(prob, progress=true, progress_steps=1)

#%% Run the analysis
print("[$(now())]: Running analysis... ")
ts, data = timeseries_analysis(sol)
burst_tstamps = ts["Bursts"][2] #Extract the bursts
println(" Completed")


#%% Extract plotting data
print("[$(now())]: Extracting data... ")
offset = 2000
tstops = burst_tstamps[2, 1]-offset:dt:burst_tstamps[2, 2]+offset
vt = sol(tstops, idxs=2).u
et = sol(tstops, idxs=9).u
it = sol(tstops, idxs=10).u
tstops = (tstops .- tstops[1]) ./ 1000

v_rng = -80.0:1.0:0.0
pVSe = parameters[VSe]
pV0e = parameters[V0e]
pVSi = parameters[VSi]
pV0i = parameters[V0i]
pρe = parameters[ρe]
pρi = parameters[ρi]
ACH_r = Φ.(v_rng, pVSe, pV0e) * pρe
GABA_r = Φ.(v_rng, pVSi, pV0i) * pρi
println(" Completed")

# Generate a release map
print("[$(now())]: Generating a neurotransmitter release map... ")
pDe = parameters[De]
pDi = parameters[Di]
nx = ny = 25
ACH_map = zeros(nx, ny, length(tstops))
GABA_map = zeros(nx, ny, length(tstops))

center_x = round(Int64, nx / 2) #Because the index starts from 0 subtract 1
center_y = round(Int64, ny / 2) #Because the indexing starts from 0
# Calculate each frames diffusion
dXe = (1.0, 1.0)
dYe = (1.0, 1.0)
dXi = (1.0, 1.0)
dYi = (0.1, 1.9)
for i = 1:length(tstops)
     if i == 1
          ACH_map[center_x+1, center_y+1, 1] = et[i] #Set intial map
          GABA_map[center_x+1, center_y+1, 1] = it[i] #Set initial map
     else
          #println(i)
          e_i = ACH_map[:, :, i-1]
          i_i = GABA_map[:, :, i-1]
          de = zeros(size(e_i))
          di = zeros(size(i_i))
          ∇(de, e_i, pDe, dX=dXe, dY=dYe)
          ∇(di, i_i, pDi, dX=dXi, dY=dYi)
          println(sum(e_i))
          ACH_map[:, :, i] += de + e_i
          GABA_map[:, :, i] += di + i_i
          ACH_map[center_x+1, center_y+1, i] = et[i] #Set intial map
          GABA_map[center_x+1, center_y+1, i] = it[i] #Set initial map
     end
end

# Set the frame stops for Images
frame_stops = LinRange(burst_tstamps[2, 1] + dt, burst_tstamps[2, 2] + 200, 4) .- burst_tstamps[2, 1] .+ (offset)
frame_stops = round.(Int64, frame_stops)
avg_ACH = sum(ACH_map, dims=3)[:, :, 1] ./ length(tstops)
avg_GABA = sum(GABA_map, dims=3)[:, :, 1] ./ length(tstops)

#Generate current maps
reload_parameters() #Might need to reload parameters
pgACh = parameters[g_ACh]
pEACh = parameters[E_ACh]
pkACh = parameters[k_ACh]
pgGABA = parameters[g_GABA]
pECl = parameters[E_Cl]
pkGABA = parameters[k_GABA]
V_CLAMP = -40.0
I_ACh = -pgACh .* ħ.(ACH_map, pkACh) .* (V_CLAMP - pEACh)
I_GABA = -pgGABA .* ħ.(GABA_map, pkGABA) .* (V_CLAMP - pECl)
I_TOTAL = I_ACh + I_GABA
avg_I_TOTAL = sum(I_TOTAL, dims=3)[:, :, 1] ./ size(I_TOTAL, 3)

#Maximal calculations
I_ACh_max = -pgACh .* 1 .* (V_CLAMP - pEACh)
I_GABA_max = -pgGABA .* 1 .* (V_CLAMP - pECl)

# Finally we need to extract some data from the X or y
dY_D = 3:2:center_y
dY_V = center_y+3:2:nx-1
dX_N = 3:2:center_x
dX_T = center_x+3:2:ny-1

#%% Lets plot
print("[$(now())]: Plotting... ")
width_inches = 7.5
height_inches = 7.5
fig3 = plt.figure("Neurotransmitter Dynamics", figsize=(width_inches, height_inches))

gs = fig3.add_gridspec(4, 2,
     width_ratios=(0.77, 0.23),
     height_ratios=(0.2, 0.26, 0.26, 0.26),
     right=0.88, left=0.12,
     top=0.95, bottom=0.08,
     wspace=0.20, hspace=0.4
)
col1_ylabel = -0.1
col2_ylabel = -0.25
im1_levels = 0.0:0.05:1.0
im2_levels = 0.0:0.05:1.0
im3_levels = -13.0:0.5:13.0
#% =====================================================Make panel A===================================================== %%#
gsA = gs[1, 1].subgridspec(ncols=1, nrows=2)

axA1 = fig3.add_subplot(gsA[1, 1])
axA1.plot(tstops, vt, c=v_color, lw=lw_standard)
ylabel("Volt. \n(mV)")
axA1.xaxis.set_visible(false) #Turn off the bottom axis
axA1.yaxis.set_label_coords(col1_ylabel, 0.5)
axA1.spines["bottom"].set_visible(false)

axA2 = fig3.add_subplot(gsA[2, 1])
axA2.plot(tstops, et, c=:green, lw=lw_standard)
axA2.plot(tstops, it, c=:red, lw=lw_standard)
axA2.set_xlabel("Time (s)")
axA2.set_ylabel("Rel. \n(mM)")
axA2.set_ylim(0.0, 2.0)
#axA2.xaxis.set_visible(false) #Turn off the bottom axis
axA2.yaxis.set_label_coords(col1_ylabel, 0.5)
axA2.legend(["ACh", "GABA"],
     ncol=2, columnspacing=0.70,
     bbox_to_anchor=(0.65, 1.1), fontsize=9.0, markerscale=0.5, handletextpad=0.3
)
#Plot all of the frame stops
axA2.plot(tstops[frame_stops], et[frame_stops], c=:blue, markersize=2.0, linewidth=0.0, marker="o")
axA2.plot(tstops[frame_stops], it[frame_stops], c=:red, markersize=2.0, linewidth=0.0, marker="o")

axAR = fig3.add_subplot(gs[1, 2])
axAR.set_xlabel("Voltage (mV)")
axAR.set_ylabel("Rel. (mM)")
axAR.plot(v_rng, ACH_r, c=:green, lw=lw_standard)
axAR.plot(v_rng, GABA_r, c=:red, lw=lw_standard)

#% ===============================================Make panel B=============================================== %%#
gsBL = gs[2, 1].subgridspec(ncols=4, nrows=1)
#We should pick 4 locations from 
for (idx, fr) in enumerate(frame_stops)
     axBi = fig3.add_subplot(gsBL[idx])
     if idx == 1
          ylabel("ACh Release")
          axBi.set_yticks([])
          axBi.yaxis.set_label_coords(col2_ylabel, 0.5)
     else
          axBi.yaxis.set_visible(false) #Turn off the bottom axis
     end
     ctr_B = axBi.contourf(ACH_map[:, :, fr], cmap="Greens", levels=im1_levels, vmin=0.0, vmax=0.50, extend="both")
     axBi.xaxis.set_visible(false) #Turn off the bottom axis
     axBi.set_aspect("equal", "box")
     axBi.spines["left"].set_visible(false)
     axBi.spines["bottom"].set_visible(false)
end

axBR = fig3.add_subplot(gs[2, 2])
axBR.set_aspect("equal", "box")
ctr_E = axBR.contourf(avg_ACH, cmap="Greens", levels=im1_levels, vmin=0.0, vmax=0.50, extend="both")
axBR.set_xticks([])
axBR.yaxis.set_visible(false) #Turn off the bottom axis
axBR.spines["bottom"].set_visible(false)
axBR.spines["left"].set_visible(false)
xlabel("Avg. ACh \n Release (mM)")
cbarE = fig3.colorbar(ctr_E, ticks=[0.0, 0.50, 1.0], aspect=5)
cbarE.ax.set_ylabel("Average Et (ACh)")


#% ===============================================Make panel C=============================================== %%#
gsCL = gs[3, 1].subgridspec(ncols=4, nrows=1)

for (idx, fr) in enumerate(frame_stops)
     axCi = fig3.add_subplot(gsCL[idx])
     ctr_I = axCi.contourf(GABA_map[:, :, fr], cmap="Reds", levels=im2_levels, vmin=0.0, vmax=0.50, extend="both")
     #axCi.xaxis.set_visible(false) #Turn off the bottom axis
     if idx == 1
          ylabel("GABA Release")
          axCi.set_yticks([])
          axCi.yaxis.set_label_coords(col2_ylabel, 0.5)
     else
          axCi.yaxis.set_visible(false) #Turn off the bottom axis
     end
     axCi.set_xticks([])
     axCi.set_aspect("equal", "box")
     xlabel("t = $(round(tstops[fr], digits = 1))")
     axCi.spines["left"].set_visible(false)
     axCi.spines["bottom"].set_visible(false)
     #Can we put a text box below each frame 

end

axCR = fig3.add_subplot(gs[3, 2])
ctr_I = axCR.contourf(avg_GABA, cmap="Reds", levels=im2_levels, vmin=0.0, vmax=0.50, extend="both")
axCR.set_xticks([])
axCR.yaxis.set_visible(false) #Turn off the bottom axis
axCR.set_aspect("equal", "box")
axCR.spines["bottom"].set_visible(false)
axCR.spines["left"].set_visible(false)
xlabel("Avg. GABA \n Release (mM)")
cbarI = fig3.colorbar(ctr_I, ticks=[0.0, 0.50, 1.0], aspect=5)
cbarI.ax.set_ylabel("Average It (GABA)")

#% ===============================================Make panel D=============================================== %%#
gsDL = gs[4, 1].subgridspec(ncols=3, nrows=3)
axD1 = fig3.add_subplot(gsDL[1, 2])
cmapYV = plt.get_cmap("Blues")
cmapYD = plt.get_cmap("Greens")
cmapXN = plt.get_cmap("Reds")
cmapXT = plt.get_cmap("Purples")
#Plot the Nasal direction
for (i, dyv) in enumerate(dY_V)
     axD1.plot(I_TOTAL[dyv, center_x, :], c=cmapYV(i / length(dY_V)), lw=lw_standard)
end
axD1.set_ylim(-15.0, 10.0)

axD1.xaxis.set_visible(false) #Turn off the bottom axis
axD1.spines["bottom"].set_visible(false)

axD2 = fig3.add_subplot(gsDL[3, 2])
#Plot the Nasal direction
for (i, dyd) in enumerate(reverse(dY_D))
     axD2.plot(tstops, I_TOTAL[dyd, center_x, :], c=cmapYD(i / length(dY_D)), lw=lw_standard)
end
axD2.set_xlabel("Time (s)")
axD2.set_ylim(-15.0, 10.0)

axDC = fig3.add_subplot(gsDL[2, 2])
axDC.plot(tstops, I_TOTAL[center_y, center_x, :], c=:black, lw=lw_standard)
axDC.set_ylim(-15.0, 10.0)
axDC.xaxis.set_visible(false) #Turn off the bottom axis
axDC.spines["bottom"].set_visible(false)
axDC.yaxis.set_visible(false)
axDC.spines["left"].set_visible(false)

axD3 = fig3.add_subplot(gsDL[2, 1])
for (i, dxn) in enumerate(reverse(dX_N))
     axD3.plot(tstops, I_TOTAL[center_y, dxn, :], c=cmapXN(i / length(dX_N)), lw=lw_standard)
end
axD3.set_xlabel("Time (s)")
axD3.set_ylim(-15.0, 10.0)
axD3.set_ylabel("Current (pA)")
#axD3.xaxis.set_visible(false) #Turn off the bottom axis
#axD3.spines["bottom"].set_visible(false)

axD4 = fig3.add_subplot(gsDL[2, 3])
ylim(-15.0, 10.0)
for (i, dxt) in enumerate(dX_T)
     axD4.plot(tstops, I_TOTAL[center_y, dxt, :], c=cmapXT(i / length(dX_T)), lw=lw_standard)
end
axD4.set_xlabel("Time (s)")

axDR = fig3.add_subplot(gs[4, 2])
ctr_i = axDR.contourf(avg_I_TOTAL, cmap="RdYlGn", levels=im3_levels, vmin=-5.0, vmax=5.0, extend="both")
axDR.plot([center_x], [center_y], linewidth=0.0, marker="o", ms=4.0)
#Plot sample points
valYD = abs.(center_y .- dY_D .- 1)
valYV = abs.(center_y .- dY_V .- 1)
axDR.scatter(fill(center_x, length(dY_D)), dY_D .- 1, s=4.0, c=valYD, cmap=cmapYD, lw=lw_standard, marker="s") #Dorsal
axDR.scatter(fill(center_x, length(dY_V)), dY_V .- 1, s=4.0, c=valYV, cmap=cmapYV, lw=lw_standard, marker="s") #Dorsal
valXN = abs.(center_x .- dX_N .- 1)
valXT = abs.(center_x .- dX_T .- 1)
axDR.scatter(dX_N .- 1, fill(center_y, length(dX_N)), s=4.0, c=valXN, cmap=cmapXN, lw=lw_standard, marker="s") #Dorsal
axDR.scatter(dX_T .- 1, fill(center_y, length(dX_T)), s=4.0, c=valXT, cmap=cmapXT, lw=lw_standard, marker="s") #Dorsal
axDR.set_facecolor("none")

axDR.set_xticks([])
axDR.yaxis.set_visible(false) #Turn off the bottom axis
axDR.set_xlim(0, nx)
axDR.set_ylim(0, ny)
axDR.set_aspect("equal", "box")
axDR.spines["bottom"].set_visible(false)
axDR.spines["left"].set_visible(false)
cbari = fig3.colorbar(ctr_i, ticks=[-10.0, 0.0, 10.0], aspect=5)
cbari.ax.set_ylabel("Induced Current (pA)")

axDR.annotate("A", (0.01, 0.95), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axDR.annotate("B", (0.01, 0.75), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axDR.annotate("C", (0.01, 0.50), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axDR.annotate("D", (0.01, 0.25), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
println(" Completed")

#%% Save the plot
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 2\Figures"
print("[$(now())]: Saving the figure 3...")
fig3.savefig("$(loc)/Figure5_Neurotransmitters.jpg")
plt.close("all")
println(" Completed")