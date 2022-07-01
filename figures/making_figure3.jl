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
import RetinalChaos: T_ODE_NT_Release
import Plots: contourf
include("figure_setup.jl")
println("Running the plotting script for figure 3")

#%% Set up modelling data
print("[$(now())]: Setting up modelling data... ")
#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
#Step 2: Import the parameters
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 15.0
pars_dict[:g_ACh] = 0.0
pars_dict[:g_GABA] = 0.0
p = pars_dict |> extract_dict
#Step 3: determine the timespan
tspan = (0.0, 300e3)
#Step 4: set up the problem
prob = ODEProblem(T_ODE_NT_Release, u0, tspan, p)
#Step 5: Solve the problem
sol = solve(prob, progress=true, progress_steps=1);
println(" Completed")

#%% Run the analysis
print("[$(now())]: Running analysis... ")
dt = 1.0; #Set the time differential
thresh = calculate_threshold(sol); #Extract the threshold
spike_tstamps = get_timestamps(sol); #Extract the spikes
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps); #Extract the bursts
println(" Completed")

#%% Extract plotting data
print("[$(now())]: Extracting data... ")
offset = 2000
t = burst_tstamps[2, 1]-offset:dt:burst_tstamps[2, 2]+offset
vt = sol(t, idxs=1)
et = sol(t, idxs=6)
it = sol(t, idxs=7)
t = (t .- t[1]) ./ 1000

v_rng = -80.0:1.0:0.0
Vse = pars_dict[:Vse]
V0e = pars_dict[:V0e]
Vsi = pars_dict[:Vsi]
V0i = pars_dict[:V0i]
ρe = pars_dict[:ρe]
ρi = pars_dict[:ρi]
ACH_r = Φ.(v_rng, Vse, V0e) * ρe
GABA_r = Φ.(v_rng, Vsi, V0i) * ρi
println(" Completed")

# Generate a release map
print("[$(now())]: Generating a neurotransmitter release map... ")
De = pars_dict[:De]
Di = pars_dict[:Di]
nx = ny = 24
ACH_map = zeros(nx, ny, length(t))
GABA_map = zeros(nx, ny, length(t))

# Calculate each frames diffusion
dXe = (1.0, 1.0)
dYe = (1.0, 1.0)
dXi = (1.0, 1.0)
dYi = (0.1, 1.9)
for i = 1:length(t)
     if i == 1
          ACH_map[round(Int64, nx / 2), round(Int64, nx / 2), 1] = et[i] #Set intial map
          GABA_map[round(Int64, nx / 2), round(Int64, nx / 2), 1] = it[i] #Set initial map
     else
          #println(i)
          e_i = ACH_map[:, :, i-1]
          i_i = GABA_map[:, :, i-1]
          de = zeros(size(e_i))
          di = zeros(size(i_i))
          ∇(de, e_i, De, dX=dXe, dY=dYe)
          ∇(di, i_i, Di, dX=dXi, dY=dYi)
          println(sum(e_i))
          ACH_map[:, :, i] += de + e_i
          GABA_map[:, :, i] += di + i_i
          ACH_map[round(Int64, nx / 2), round(Int64, nx / 2), i] = et[i] #Set intial map
          GABA_map[round(Int64, nx / 2), round(Int64, nx / 2), i] = it[i] #Set initial map
     end
end

# Set the frame stops for Images
frame_stops = LinRange(burst_tstamps[2, 1] + dt, burst_tstamps[2, 2] + 200, 4) .- burst_tstamps[2, 1] .+ (offset)
frame_stops = round.(Int64, frame_stops)
avg_ACH = sum(ACH_map, dims=3)[:, :, 1] ./ length(t)
avg_GABA = sum(GABA_map, dims=3)[:, :, 1] ./ length(t)

#Generate current maps
pars_dict = read_JSON("params\\params.json") #Might need to reload parameters
gACh = pars_dict[:g_ACh]
EACh = pars_dict[:E_ACh]
kACh = pars_dict[:k_ACh]
gGABA = pars_dict[:g_GABA]
EGABA = pars_dict[:E_GABA]
kGABA = pars_dict[:k_GABA]
V_CLAMP = -50.0
I_ACh = -gACh .* ħ.(ACH_map, kACh) .* (V_CLAMP - EACh)
I_GABA = -gGABA .* ħ.(GABA_map, kACh) .* (V_CLAMP - EGABA)
I_TOTAL = I_ACh + I_GABA
avg_I_TOTAL = sum(I_TOTAL, dims=3)[:, :, 1] ./ size(I_TOTAL, 3)

# Finally we need to extract some data from the X or y
center_x = round(Int64, nx / 2)
center_y = round(Int64, ny / 2)
dX_N = round.(Int64, LinRange(1, center_x-1, 8))
dX_T = round.(Int64, LinRange(center_x, nx, 8))
dY_V = round.(Int64, LinRange(1, center_y-1, 8))
dY_D = round.(Int64, LinRange(center_y, ny, 8))

#%% Lets plot
print("[$(now())]: Plotting... ")
width_inches = 16.0
height_inches = 12.0
fig3 = plt.figure("Neurotransmitter Dynamics", figsize=(width_inches, height_inches))

gs = fig3.add_gridspec(4, 2,
     width_ratios=(0.77, 0.23),
     height_ratios=(0.2, 0.26, 0.26, 0.26),
     right=0.95, left=0.07,
     top=0.95, bottom=0.08,
     wspace=0.10, hspace=0.4
)
col1_ylabel = -0.05
col2_ylabel = -0.17
#% =====================================================Make panel A===================================================== %%#
gsA = gs[1, 1].subgridspec(ncols=1, nrows=2)

axA1 = fig3.add_subplot(gsA[1, 1])
axA1.plot(t, vt, c=v_color, lw=3.0)
ylabel("Vt (mV)")
axA1.xaxis.set_visible(false) #Turn off the bottom axis
axA1.yaxis.set_label_coords(col1_ylabel, 0.5)

axA2 = fig3.add_subplot(gsA[2, 1])
ylim(0.0, 5.0)
axA2.plot(t, et, c=e_color, lw=3.0)
ylabel("Et")
axA2.xaxis.set_visible(false) #Turn off the bottom axis
axA2.yaxis.set_label_coords(col1_ylabel, 0.5)

axA3 = fig3.add_subplot(gsA[2, 1])
ylim(0.0, 5.0)
axA3.plot(t, it, c=:red, lw=3.0)
xlabel("Time (s)")
ylabel("It")
axA3.yaxis.set_label_coords(col1_ylabel, 0.5)
#Plot all of the frame stops
axA3.plot(t[frame_stops], et[frame_stops], c=:blue, linewidth=0.0, marker="o")
axA3.plot(t[frame_stops], it[frame_stops], c=:red, linewidth=0.0, marker="o")

axAR = fig3.add_subplot(gs[1, 2])
axAR.plot(v_rng, ACH_r, c=e_color, lw=3.0)
axAR.plot(v_rng, GABA_r, c=:red, lw=3.0)



#% ===============================================Make panel B=============================================== %%#
gsBL = gs[2, 1].subgridspec(ncols=4, nrows=1)
#We should pick 4 locations from 
for (idx, fr) in enumerate(frame_stops)
     axBi = fig3.add_subplot(gsBL[idx])
     ctr_B = axBi.contourf(ACH_map[:, :, fr], cmap="inferno", levels=0.0:0.05:0.5, extend="both")
     axBi.xaxis.set_visible(false) #Turn off the bottom axis
     axBi.yaxis.set_visible(false) #Turn off the bottom axis
     axBi.set_aspect("equal", "box")
end

axBR = fig3.add_subplot(gs[2, 2])
axBR.set_aspect("equal", "box")
ctr_E = axBR.contourf(avg_ACH, cmap="inferno", levels=0.0:0.05:0.5, extend="both")
axBR.set_xticks([])
axBR.yaxis.set_visible(false) #Turn off the bottom axis
xlabel("Avg. ACh Release")
cbarE = fig3.colorbar(ctr_E, ticks=[0.0, 0.25, 0.5])
cbarE.ax.set_ylabel("Average Et (ACh)")


#% ===============================================Make panel C=============================================== %%#
gsCL = gs[3, 1].subgridspec(ncols=4, nrows=1)

for (idx, fr) in enumerate(frame_stops)
     axCi = fig3.add_subplot(gsCL[idx])
     ctr_I = axCi.contourf(GABA_map[:, :, fr], cmap="inferno", levels=0.0:0.05:0.5, extend="both")
     #axCi.xaxis.set_visible(false) #Turn off the bottom axis
     axCi.set_xticks([])
     axCi.yaxis.set_visible(false) #Turn off the bottom axis
     axCi.set_aspect("equal", "box")
     xlabel("t = $(round(t[fr], digits = 1))")
     #Can we put a text box below each frame 

end

axCR = fig3.add_subplot(gs[3, 2])
ctr_I = axCR.contourf(avg_GABA, cmap="inferno", levels=0.0:0.05:0.5, extend="both")
axCR.set_xticks([])
axCR.yaxis.set_visible(false) #Turn off the bottom axis
axCR.set_aspect("equal", "box")
xlabel("Avg. GABA Release")
cbarI = fig3.colorbar(ctr_I, ticks=[0.0, 0.25, 0.5])
cbarI.ax.set_ylabel("Average It (GABA)")
#% ===============================================Make panel D=============================================== %%#
gsDL = gs[4, 1].subgridspec(ncols=2, nrows=2)
axD1 = fig3.add_subplot(gsDL[1, 1])
ylim(-10.0, 10.0)
#Plot the Nasal direction
for dxn in dX_N
     axD1.plot(I_TOTAL[dxn, center_y, :])
end
axD1.xaxis.set_visible(false) #Turn off the bottom axis
axD1.spines["bottom"].set_visible(false)

axD2 = fig3.add_subplot(gsDL[2, 1])
ylim(-10.0, 10.0)
#Plot the Nasal direction
for dxt in dX_T
     axD2.plot(I_TOTAL[dxt, center_y, :])
end

axD3 = fig3.add_subplot(gsDL[1, 2])
ylim(-10.0, 10.0)
for dyd in dY_D
     axD3.plot(I_TOTAL[center_x, dyd, :])
end
axD3.xaxis.set_visible(false) #Turn off the bottom axis
axD3.spines["bottom"].set_visible(false)

axD4 = fig3.add_subplot(gsDL[2, 2])
ylim(-10.0, 10.0)
for dyv in dY_V
     axD4.plot(I_TOTAL[center_x, dyv, :])
end

axDR = fig3.add_subplot(gs[4, 2])
ctr_i = axDR.contourf(avg_I_TOTAL, cmap="RdYlGn", levels=-5.0:0.5:5.0, extend="both")
axDR.set_xticks([])
axDR.yaxis.set_visible(false) #Turn off the bottom axis
axDR.set_aspect("equal", "box")
axDR.spines["bottom"].set_visible(false)
axDR.spines["left"].set_visible(false)
xlabel("Induced Current")
cbari = fig3.colorbar(ctr_i, ticks=[-10.0, 0.0, 10.0])
cbari.ax.set_ylabel("Current (pA)")
#Lets put text labeling nasal temporal 
axDR.text(0.0, 10.0, "N")
axDR.text(20.0, 10.0, "T")
axDR.text(10.0, 0.0, "D")
axDR.text(10.0, 20.0, "V")
println(" Completed")






#%% Save the plot
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2021 A Computational Model - Sci. Rep\Figures"
print("[$(now())]: Saving the figure 2...")
fig3.savefig("$(loc)/figure2_BiophysicalProperties.png")
plt.close("All")
println(" Completed")

#%% Generate a diffusion animation
anim = @animate for i = 1:100:length(t)
     print("[$(now())]: Animating simulation...")
     println(i)
     e_frame = ACH_map[:, :, i]
     i_frame = GABA_map[:, :, i]
     I_frame = I_TOTAL[:, :, i]
     #contourf(frame_e, e_frame, clims=(0.0, 1.0))
     contour_e = contourf(e_frame,
          ratio=:equal, grid=false,
          ylabel="t = $(t[i])",
          xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
          c=:curl, clims=(0.0, 1.0), levels=0.0:0.1:1.0
     )
     contour_i = contourf(i_frame,
          ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
          c=:curl, clims=(0.0, 1.0), levels=0.0:0.1:1.0
     )

     contour_I = contourf(I_frame,
          ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
          c=:RdBu, clims=(-10.0, 10.0), levels=-10.0:0.5:10.0
     )
     plot(contour_e, contour_i, contour_I, layout=3)
end

gif(anim, "$(loc)/supplemental_animation.gif")