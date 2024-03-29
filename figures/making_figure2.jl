include("figure_setup.jl")
using RetinalChaos
import RetinalChaos: get_timestamps
plt.pygui(true) #Make the GUI external to vscode
println("Running the plotting script for figure 1")
#%% Open data
#Step 1: Import the model, initial conditions, and parameters
import RetinalChaos.ODEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters

#Step 2: Adjust parameters
parameters[I_app] = 15.0
parameters[ρi] = 0.0
parameters[ρe] = 0.0
parameters[g_Na] = 0.0 #Don't include voltage gated sodium channels

#Step 3: determine the timespan
tmin = 0.0
tmax = 120e3

#Step 4: set up the problem
prob = ODEProblem(ODEModel, u0, (tmin, tmax), parameters)

#Step 5: Solve the problem
sol = solve(prob, progress=true, progress_steps=1);
println(" Completed")
ODEModel.states

fig, ax = plt.subplots(1)
ax.plot(sol.t, sol(sol.t, idxs = v))
fig
# Run the analysis 
print("[$(now())]: Running analysis... ")
ts, data = RetinalChaos.timeseries_analysis(sol)
tSpike, vnSpike = extract_spike_trace(ts, data, cell_n = 2, idx = 1, spike_dur=200)
ts["Spikes"]

#Upsample the spike time
dtSpike = 0.01 #Set the time differential
tSpike = tSpike[1]:dtSpike:tSpike[end]
vSpike = sol(tSpike, idxs=2)
nSpike = sol(tSpike, idxs=3)
#tSpike = tSpike .- tSpike[1]

#Spikes arrows
deltaArrow = 1.0
tSpikeArrow = LinRange(tSpike[1], tSpike[end], 25)
vSpikeArrow = sol(tSpikeArrow, idxs=[2]) |> Array
nSpikeArrow = sol(tSpikeArrow, idxs=[3]) |> Array
dvSpikeArrow = sol(tSpikeArrow .+ deltaArrow, idxs=[1]) .- vSpikeArrow
dnSpikeArrow = sol(tSpikeArrow .+ deltaArrow, idxs=[2]) .- nSpikeArrow

#extract the bursts
beginBurst = ts["Bursts"][2][2]
dtBurst = 1.0
tBurst = (beginBurst-1000):dtBurst:(beginBurst+2000)
vBurst = sol(tBurst, idxs=2)
cBurst = sol(tBurst, idxs=6)
#tBurst = (tBurst .- tBurst[1]) ./ 1000 #Offset the time range

#Bursts
tBurstArrow = tBurst[1]:250:tBurst[end]
vBurstArrow = sol(tBurstArrow, idxs=[2]) |> Array
cBurstArrow = sol(tBurstArrow, idxs=[6]) |> Array
dvBurstArrow = sol(tBurstArrow .+ deltaArrow, idxs=[1]) .- vBurstArrow
dcBurstArrow = sol(tBurstArrow .+ deltaArrow, idxs=[3]) .- cBurstArrow

#extract the IBI
C_dt = 1.0
ts["Bursts"][2][2,2]
tIBI = (ts["Bursts"][2][2, 2]-1000):C_dt:(ts["Bursts"][2][3, 2]+10000)
vIBI = sol(tIBI, idxs=2)
cIBI = sol(tIBI, idxs=6)
aIBI = sol(tIBI, idxs=7)
bIBI = sol(tIBI, idxs=8)

#IBIs
tIBIArrow = LinRange(tIBI[1], tIBI[end], 100)
vIBIArrow = sol(tIBIArrow, idxs=[2]) |> Array
cIBIArrow = sol(tIBIArrow, idxs=[6]) |> Array
aIBIArrow = sol(tIBIArrow, idxs=[7]) |> Array
bIBIArrow = sol(tIBIArrow, idxs=[8]) |> Array

dvIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[2]) .- vIBIArrow
dcIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[6]) .- cIBIArrow
daIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[7]) .- cIBIArrow
dbIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[8]) .- cIBIArrow
println(" Completed")

#plt.clf() #While drawing you can use this to clear the figure 

#%% Start the plotting of the figure
print("[$(now())]: Plotting figure 1...")
width_inches = 7.5
height_inches = 7.5
fig1 = plt.figure("Model Basics", figsize=(width_inches, height_inches))

#% Make a plot in PyPlot
gs = fig1.add_gridspec(3, 2,
     width_ratios=(0.75, 0.25),
     height_ratios=(0.30, 0.30, 0.40),
     right=0.95, left=0.14,
     top=0.94, bottom=0.08,
     wspace=0.4, hspace=0.40
)

col1_ylabel = -0.1
col2_ylabel = -0.25
#% =====================================================Make panel A===================================================== %%#
vlims = (-55.0, 5.0)
nlims = (-0.1, 1.1)
clims = (0.0, 0.5)
alims = (-0.1, 0.6)
blims = (-0.1, 0.6)

#Do the plotting
gsAL = gs[1, 1].subgridspec(ncols=1, nrows=2)
axA1 = fig1.add_subplot(gsAL[1])
ylim(vlims)
axA1.plot(tSpike .- tSpike[1], vSpike, c=v_color, lw=lw_standard)
ylabel("Voltage \n (mV)")
axA1.yaxis.set_label_coords(col1_ylabel, 0.5)
axA1.xaxis.set_visible(false) #Turn off the bottom axis
axA1.spines["bottom"].set_visible(false)
axA1.yaxis.set_major_locator(MultipleLocator(50.0))
axA1.yaxis.set_minor_locator(MultipleLocator(25.0))

axA2 = fig1.add_subplot(gsAL[2])
ylim(nlims)
axA2.plot(tSpike .- tSpike[1], nSpike, c=n_color, lw=lw_standard)
ylabel("K-Chan. \n open prob.")
xlabel("Time (ms)")
axA2.yaxis.set_label_coords(col1_ylabel, 0.5)
axA2.yaxis.set_major_locator(MultipleLocator(1.00))
axA2.yaxis.set_minor_locator(MultipleLocator(0.50))

axAR = fig1.add_subplot(gs[1, 2])
xlim(vlims)
ylim(nlims)
axAR.plot(vSpike, nSpike, c=:black, lw=lw_standard)
add_direction(axAR, vSpikeArrow, nSpikeArrow, dvSpikeArrow, dnSpikeArrow)
xlabel("Voltage (mV)")
ylabel("K-Chan. \n open prob.")
axAR.xaxis.set_major_locator(MultipleLocator(25.0))
axAR.xaxis.set_minor_locator(MultipleLocator(12.5))
axAR.yaxis.set_major_locator(MultipleLocator(0.50))
axAR.yaxis.set_minor_locator(MultipleLocator(0.25))
axAR.yaxis.set_label_coords(col2_ylabel, 0.5)
#add arrows here?
#% ===============================================Make panel B=============================================== %%#
vlims = (-75.0, 5.0)
gsBL = gs[2, 1].subgridspec(ncols=1, nrows=2)
axB1 = fig1.add_subplot(gsBL[1])
ylim(vlims)
axB1.plot((tBurst .- tBurst[1]) ./ 1000, vBurst, c=v_color, lw=lw_standard)
axB1.xaxis.set_visible(false) #Turn off the bottom axis
ylabel("Voltage \n (mV)")
axB1.yaxis.set_label_coords(col1_ylabel, 0.5)
axB1.spines["bottom"].set_visible(false)
axB1.yaxis.set_major_locator(MultipleLocator(30.0))
axB1.yaxis.set_minor_locator(MultipleLocator(15.0))

axB2 = fig1.add_subplot(gsBL[2])
ylim(0.0, 0.41)
axB2.plot((tBurst .- tBurst[1]) ./ 1000, cBurst, c=c_color, lw=lw_standard)
xlabel("Time (s)")
ylabel("[Ca2+] \n (mM)")
axB2.yaxis.set_label_coords(col1_ylabel, 0.5)
axB2.yaxis.set_major_locator(MultipleLocator(0.20))
axB2.yaxis.set_minor_locator(MultipleLocator(0.1))

axBR = fig1.add_subplot(gs[2, 2])
xlim(-75.0, 5.0)
ylim(0.0, 0.41)
add_direction(axBR, vBurstArrow, cBurstArrow, dvBurstArrow, dcBurstArrow, color=c_color)
axBR.plot(vBurst, cBurst, c=c_color, lw=2.0)
xlabel("Voltage (mV)")
ylabel("[Ca2+] (mM)")
axBR.yaxis.set_label_coords(col2_ylabel, 0.5)
axBR.xaxis.set_major_locator(MultipleLocator(30.0))
axBR.xaxis.set_minor_locator(MultipleLocator(15.0))
axBR.yaxis.set_major_locator(MultipleLocator(0.2))
axBR.yaxis.set_minor_locator(MultipleLocator(0.1))

#% ===================================================Figure Panel C=================================================== %%#
gsCL = gs[3, 1].subgridspec(ncols=1, nrows=4)
axCL1 = fig1.add_subplot(gsCL[1])
ylim(-0.1, 0.41)
axCL1.plot((tIBI .- tIBI[1]) ./ 1000, cIBI, c=c_color, lw=lw_standard)
ylabel("[Ca2+] \n (mM)")
axCL1.xaxis.set_visible(false) #Turn off the bottom axis
axCL1.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL1.spines["bottom"].set_visible(false)
axCL1.yaxis.set_major_locator(MultipleLocator(0.4))
axCL1.yaxis.set_minor_locator(MultipleLocator(0.2))

axCL2 = fig1.add_subplot(gsCL[2])
ylim(-0.76, 0.1)
axCL2.plot((tIBI .- tIBI[1]) ./ 1000, -aIBI, c=a_color, lw=lw_standard)
ylabel("cAMP \n (At)")
axCL2.xaxis.set_visible(false) #Turn off the bottom axis
axCL2.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL2.spines["bottom"].set_visible(false)
axCL2.yaxis.set_major_locator(MultipleLocator(0.60))
axCL2.yaxis.set_minor_locator(MultipleLocator(0.30))

axCL3 = fig1.add_subplot(gsCL[3])
ylim(-0.1, 0.76)
axCL3.plot((tIBI .- tIBI[1]) ./ 1000, bIBI, c=b_color, lw=lw_standard)
ylabel("TREK \n (Bt)")
axCL3.xaxis.set_visible(false) #Turn off the bottom axis
axCL3.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL3.spines["bottom"].set_visible(false)
axCL3.yaxis.set_major_locator(MultipleLocator(0.60))
axCL3.yaxis.set_minor_locator(MultipleLocator(0.30))

axCL4 = fig1.add_subplot(gsCL[4])
ylim(-90.0, 5.0)
axCL4.plot((tIBI .- tIBI[1]) ./ 1000, vIBI, c=v_color, lw=lw_standard)
xlabel("Time (s)")
ylabel("Volt. \n (mV)")
axCL4.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL4.yaxis.set_major_locator(MultipleLocator(60.0))
axCL4.yaxis.set_minor_locator(MultipleLocator(30.0))

gsCR = gs[3, 2].subgridspec(ncols=1, nrows=5, height_ratios=[0.2, 0.2, 0.2, 0.2, 0.2])
axCR1 = fig1.add_subplot(gsCR[1])
xlim(0.0, 0.41)
ylim(-0.1, 0.761)
axCR1.plot(cIBI, aIBI, c=a_color, lw=lw_standard)
#add_direction(axCR1, cIBIArrow, aIBIArrow, dcIBIArrow, daIBIArrow, color=a_color)
xlabel("[Ca2+](mM)")
ylabel("cAMP \n decay (At)")
axCR1.yaxis.set_label_coords(col2_ylabel, 0.5)
axCR1.xaxis.set_major_locator(MultipleLocator(0.2))
axCR1.xaxis.set_minor_locator(MultipleLocator(0.1))
axCR1.yaxis.set_major_locator(MultipleLocator(0.60))
axCR1.yaxis.set_minor_locator(MultipleLocator(0.30))


axCR2 = fig1.add_subplot(gsCR[3])
xlim(-0.1, 0.65)
ylim(-0.1, 0.65)
axCR2.plot(aIBI, bIBI, c=b_color, lw=lw_standard)
#add_direction(axCR2, aIBIArrow, bIBIArrow, daIBIArrow, dbIBIArrow, color=b_color)
xlabel("cAMP decay (At)")
ylabel("TREK \n (Bt)")
axCR2.yaxis.set_label_coords(col2_ylabel, 0.5)
axCR2.xaxis.set_major_locator(MultipleLocator(0.3))
axCR2.xaxis.set_minor_locator(MultipleLocator(0.15))
axCR2.yaxis.set_major_locator(MultipleLocator(0.6))
axCR2.yaxis.set_minor_locator(MultipleLocator(0.3))

axCR3 = fig1.add_subplot(gsCR[5])
xlim(-0.1, 0.61)
ylim(-90.0, 5.0)
axCR3.plot(bIBI, vIBI, c=v_color, lw=lw_standard)
#add_direction(axCR3, bIBIArrow, vIBIArrow, dbIBIArrow, dvIBIArrow, color=v_color)
xlabel("TREK (Bt)")
ylabel("Voltage \n (mV)")
axCR3.yaxis.set_label_coords(col2_ylabel, 0.5)
axCR3.xaxis.set_major_locator(MultipleLocator(0.3))
axCR3.xaxis.set_minor_locator(MultipleLocator(0.15))
axCR3.yaxis.set_major_locator(MultipleLocator(60.0))
axCR3.yaxis.set_minor_locator(MultipleLocator(30.0))

axA1.annotate("A", (0.01, 0.94), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axA1.annotate("B", (0.67, 0.94), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

axA1.annotate("C", (0.01, 0.66), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axA1.annotate("D", (0.67, 0.66), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

axA1.annotate("E", (0.01, 0.35), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axA1.annotate("F", (0.67, 0.35), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

println(" Complete")
#%% Save the figure
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
print("[$(now())]: Saving the figure 1...")
fig1
fig1.savefig("$(loc)/Figure2_ModelVariables.jpg")
plt.close("all")
println(" Completed")