using RetinalChaos
import RetinalChaos: get_timestamps
include("figure_setup.jl")

println("Running the plotting script for figure 1")
#%% Open data
#Step 1: Import the initial conditions
print("[$(now())]: Setting up modelling data... ")
conds_dict = read_JSON("params\\conds.json")
conds_dict[:v] = -65.0
u0 = conds_dict |> extract_dict
#Step 2: Import the parameters
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 15.0
pars_dict[:ρi] = 0.0
pars_dict[:ρe] = 0.0
p = pars_dict |> extract_dict
#Step 3: determine the timespan
tspan = (0.0, 120e3)
#Step 4: set up the problem
prob = ODEProblem(T_ODE, u0, tspan, p)
#Step 5: Solve the problem
sol = solve(prob, progress=true, progress_steps=1);
println(" Completed")
#plot(sol)

# Run the analysis 
print("[$(now())]: Running analysis... ")
dt = 0.01 #Set the time differential
thresh = calculate_threshold(sol) #Extract the threshold
spike_tstamps = get_timestamps(sol)
spike_durs, isi = extract_interval(spike_tstamps, max_duration=100, max_interval=100)
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
burst_durs, ibi = extract_interval(burst_tstamps)
burst_tstamps
#t_series = tspan[1]:dt:tspan[end]
#vt = sol(t_series, idxs=[1])
burst_lims = burst_tstamps[2, :]#lets set the limits for area of interest
t_rng = burst_lims[1]:dt:burst_lims[2] #Set up the plotting range
println(" Completed")

# Load the data into the correct format
print("[$(now())]: Extracting data... ")
A_dx = 20 #the tick interval is 20ms
A_dt = 0.01
A_trng = (burst_lims[1]:A_dt:burst_lims[1]+200) #Set the range of points to plot
tSpike = (A_trng .- A_trng[1]) #Create the formatted time span
vSpike = sol(A_trng, idxs=1)
nSpike = sol(A_trng, idxs=2)
maximum(vSpike)
minimum(vSpike)
B_dt = 1.0
B_trng = (burst_lims[1]-1500):dt:(burst_lims[2]+1500)
tBurst = (B_trng .- B_trng[1]) ./ 1000 #Offset the time range
vBurst = sol(B_trng, idxs=1)
cBurst = sol(B_trng, idxs=3)
import RetinalChaos.M_INF
f(v) = pars_dict[:δ] * (-pars_dict[:g_Ca] * M_INF(v, pars_dict[:V1], pars_dict[:V2]) * (v - pars_dict[:E_Ca]))
vrng = -80.0:1.0:0.0
crng = f.(vrng)

C_dt = 10.0
C_trng = (burst_tstamps[2, 2]-1000):C_dt:(burst_tstamps[3, 2]+10000)

tIBI = (C_trng .- C_trng[1]) ./ 1000
vIBI = sol(C_trng, idxs=1)
cIBI = sol(C_trng, idxs=3)
aIBI = sol(C_trng, idxs=4)
bIBI = sol(C_trng, idxs=5)
println(" Completed")

#%% Prepare the arrows

#Spikes
deltaArrow = 1.0
tSpikeArrow = LinRange(A_trng[1], A_trng[end], 25)
vSpikeArrow = sol(tSpikeArrow, idxs=[1]) |> Array
nSpikeArrow = sol(tSpikeArrow, idxs=[2]) |> Array
dvSpikeArrow = sol(tSpikeArrow .+ deltaArrow, idxs=[1]) .- vSpikeArrow
dnSpikeArrow = sol(tSpikeArrow .+ deltaArrow, idxs=[2]) .- nSpikeArrow


#Bursts
#tBurstArrow = LinRange(B_trng[1], B_trng[end], 25)
tBurstArrow = B_trng[1]:250:B_trng[end]

vBurstArrow = sol(tBurstArrow, idxs=[1]) |> Array
cBurstArrow = sol(tBurstArrow, idxs=[3]) |> Array
dvBurstArrow = sol(tBurstArrow .+ deltaArrow, idxs=[1]) .- vBurstArrow
dcBurstArrow = sol(tBurstArrow .+ deltaArrow, idxs=[3]) .- cBurstArrow

#IBIs
tIBIArrow = LinRange(C_trng[1], C_trng[end], 100)
vIBIArrow = sol(tIBIArrow, idxs=[1]) |> Array
cIBIArrow = sol(tIBIArrow, idxs=[3]) |> Array
aIBIArrow = sol(tIBIArrow, idxs=[4]) |> Array
bIBIArrow = sol(tIBIArrow, idxs=[5]) |> Array

dvIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[1]) .- vIBIArrow
dcIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[3]) .- cIBIArrow
daIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[4]) .- cIBIArrow
dbIBIArrow = sol(tIBIArrow .+ deltaArrow, idxs=[5]) .- cIBIArrow

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
clims = (0.0, 0.6)
alims = (-0.1, 1.1)
blims = (-0.1, 1.1)

#Do the plotting
gsAL = gs[1, 1].subgridspec(ncols=1, nrows=2)
axA1 = fig1.add_subplot(gsAL[1])
ylim(vlims)
axA1.plot(tSpike, vSpike, c=v_color, lw=lw_standard)
ylabel("Voltage \n (mV)")
axA1.yaxis.set_label_coords(col1_ylabel, 0.5)
axA1.xaxis.set_visible(false) #Turn off the bottom axis
axA1.spines["bottom"].set_visible(false)
axA1.yaxis.set_major_locator(MultipleLocator(50.0))
axA1.yaxis.set_minor_locator(MultipleLocator(25.0))

axA2 = fig1.add_subplot(gsAL[2])
ylim(nlims)
axA2.plot(tSpike, nSpike, c=n_color, lw=lw_standard)
ylabel("K-Chan. \n open ratio")
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
ylabel("K-Chan. prob.")
axAR.xaxis.set_major_locator(MultipleLocator(25.0))
axAR.xaxis.set_minor_locator(MultipleLocator(12.5))
axAR.yaxis.set_major_locator(MultipleLocator(0.50))
axAR.yaxis.set_minor_locator(MultipleLocator(0.25))
axAR.yaxis.set_label_coords(col2_ylabel, 0.5)
#add arrows here?
#% ===============================================Make panel B=============================================== %%#
vlims = (-105.0, 5.0)

gsBL = gs[2, 1].subgridspec(ncols=1, nrows=2)
axB1 = fig1.add_subplot(gsBL[1])
ylim(vlims)
axB1.plot(tBurst, vBurst, c=v_color, lw=lw_standard)
axB1.xaxis.set_visible(false) #Turn off the bottom axis
ylabel("Voltage \n (mV)")
axB1.yaxis.set_label_coords(col1_ylabel, 0.5)
axB1.spines["bottom"].set_visible(false)
axB1.yaxis.set_major_locator(MultipleLocator(50.0))
axB1.yaxis.set_minor_locator(MultipleLocator(25.0))

axB2 = fig1.add_subplot(gsBL[2])
ylim(clims)
axB2.plot(tBurst, cBurst, c=c_color, lw=lw_standard)
xlabel("Time (s)")
ylabel("[Ca2+] \n (mM)")
axB2.yaxis.set_label_coords(col1_ylabel, 0.5)
axB2.yaxis.set_major_locator(MultipleLocator(0.25))
axB2.yaxis.set_minor_locator(MultipleLocator(0.50))

axBR = fig1.add_subplot(gs[2, 2])
xlim(vlims)
ylim(clims)
size(cBurst)
add_direction(axBR, vBurstArrow, cBurstArrow, dvBurstArrow, dcBurstArrow, color=c_color)
axBR.plot(vBurst, cBurst, c=c_color, lw=2.0)
xlabel("Voltage (mV)")
ylabel("[Ca2+] (mM)")
axBR.yaxis.set_label_coords(col2_ylabel, 0.5)
axBR.xaxis.set_major_locator(MultipleLocator(50.0))
axBR.xaxis.set_minor_locator(MultipleLocator(25.5))
axBR.yaxis.set_major_locator(MultipleLocator(0.3))
axBR.yaxis.set_minor_locator(MultipleLocator(0.15))

#% ===================================================Figure Panel C=================================================== %%#
gsCL = gs[3, 1].subgridspec(ncols=1, nrows=4)
axCL1 = fig1.add_subplot(gsCL[1])
ylim(-0.1, 0.61)
axCL1.plot(tIBI, cIBI, c=c_color, lw=lw_standard)
ylabel("[Ca2+] \n (mM)")
axCL1.xaxis.set_visible(false) #Turn off the bottom axis
axCL1.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL1.spines["bottom"].set_visible(false)
axCL1.yaxis.set_major_locator(MultipleLocator(0.6))
axCL1.yaxis.set_minor_locator(MultipleLocator(0.3))

axCL2 = fig1.add_subplot(gsCL[2])
ylim(-1.1, 0.1)
axCL2.plot(tIBI, -aIBI, c=a_color, lw=lw_standard)
ylabel("cAMP \n (At)")
axCL2.xaxis.set_visible(false) #Turn off the bottom axis
axCL2.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL2.spines["bottom"].set_visible(false)
axCL2.yaxis.set_major_locator(MultipleLocator(1.0))
axCL2.yaxis.set_minor_locator(MultipleLocator(0.5))

axCL3 = fig1.add_subplot(gsCL[3])
ylim(-0.1, 1.1)
axCL3.plot(tIBI, bIBI, c=b_color, lw=lw_standard)
ylabel("TREK \n (Bt)")
axCL3.xaxis.set_visible(false) #Turn off the bottom axis
axCL3.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL3.spines["bottom"].set_visible(false)
axCL3.yaxis.set_major_locator(MultipleLocator(1.0))
axCL3.yaxis.set_minor_locator(MultipleLocator(0.5))

axCL4 = fig1.add_subplot(gsCL[4])
ylim(-90.0, 5.0)
axCL4.plot(tIBI, vIBI, c=v_color, lw=lw_standard)
xlabel("Time (s)")
ylabel("Voltage \n (mV)")
axCL4.yaxis.set_label_coords(col1_ylabel, 0.5)
axCL4.yaxis.set_major_locator(MultipleLocator(60.0))
axCL4.yaxis.set_minor_locator(MultipleLocator(30.0))

gsCR = gs[3, 2].subgridspec(ncols=1, nrows=5, height_ratios=[0.2, 0.2, 0.2, 0.2, 0.2])
axCR1 = fig1.add_subplot(gsCR[1])
xlim(-0.1, 0.61)
ylim(-0.1, 1.1)
axCR1.plot(cIBI, aIBI, c=a_color, lw=lw_standard)
add_direction(axCR1, cIBIArrow, aIBIArrow, dcIBIArrow, daIBIArrow, color=a_color)
xlabel("[Ca2+](mM)")
ylabel("cAMP \n decay (At)")
axCR1.yaxis.set_label_coords(col2_ylabel, 0.5)
axCR1.xaxis.set_major_locator(MultipleLocator(0.5))
axCR1.xaxis.set_minor_locator(MultipleLocator(0.25))
axCR1.yaxis.set_major_locator(MultipleLocator(1.0))
axCR1.yaxis.set_minor_locator(MultipleLocator(0.5))

axCR2 = fig1.add_subplot(gsCR[3])
xlim(-0.1, 1.1)
ylim(-0.1, 1.1)
axCR2.plot(aIBI, bIBI, c=b_color, lw=lw_standard)
add_direction(axCR2, aIBIArrow, bIBIArrow, daIBIArrow, dbIBIArrow, color=b_color)
xlabel("cAMP decay (At)")
ylabel("TREK \n (Bt)")
axCR2.yaxis.set_label_coords(col2_ylabel, 0.5)
axCR2.xaxis.set_major_locator(MultipleLocator(1.0))
axCR2.xaxis.set_minor_locator(MultipleLocator(0.5))
axCR2.yaxis.set_major_locator(MultipleLocator(1.0))
axCR2.yaxis.set_minor_locator(MultipleLocator(0.5))

axCR3 = fig1.add_subplot(gsCR[5])
xlim(-0.1, 1.1)
ylim(-90.0, 5.0)
axCR3.plot(bIBI, vIBI, c=v_color, lw=lw_standard)
add_direction(axCR3, bIBIArrow, vIBIArrow, dbIBIArrow, dvIBIArrow, color=v_color)
xlabel("TREK (Bt)")
ylabel("Volt. \n (mV)")
axCR3.yaxis.set_label_coords(col2_ylabel, 0.5)
axCR3.xaxis.set_major_locator(MultipleLocator(1.0))
axCR3.xaxis.set_minor_locator(MultipleLocator(0.5))
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
fig1.savefig("$(loc)/figure1_ModelVariables.png")
plt.close("all")
println(" Completed")