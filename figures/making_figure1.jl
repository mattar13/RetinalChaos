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
plot(sol)
#%% Run the analysis 
print("[$(now())]: Running analysis... ")
dt = 0.01 #Set the time differential
thresh = calculate_threshold(sol) #Extract the threshold
spike_tstamps = get_timestamps(sol)
spike_durs, isi = extract_interval(spike_tstamps, max_duration=100, max_interval=100)
burst_tstamps, SPB = max_interval_algorithim(spike_tstamps)
burst_durs, ibi = extract_interval(burst_tstamps)
burst_tstamps
t_series = tspan[1]:dt:tspan[end]
vt = sol(t_series, idxs=[1])
burst_lims = burst_tstamps[2, :]#lets set the limits for area of interest
t_rng = burst_lims[1]:dt:burst_lims[2] #Set up the plotting range
println(" Completed")

#%% Load the data into the correct format
print("[$(now())]: Extracting data... ")
A_dx = 20 #the tick interval is 20ms
A_dt = 0.01
A_trng = (burst_lims[1]:A_dt:burst_lims[1]+200) #Set the range of points to plot
tSpike = (A_trng .- A_trng[1]) #Create the formatted time span
vSpike = sol(A_trng, idxs=1)
nSpike = sol(A_trng, idxs=2)

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

#plt.clf() #While drawing you can use this to clear the figure 
#%% Start the plotting of the figure
print("[$(now())]: Plotting figure 1...")
width_inches = 16.0
height_inches = 10.0
fig1 = plt.figure("Model Basics", figsize=(width_inches, height_inches))
fig1.text(0.0, 0.0, "A", ha="center", va="center", fontsize=12.0)

#% Make a plot in PyPlot
gs = fig1.add_gridspec(3, 2,
     width_ratios=(0.80, 0.20),
     height_ratios=(0.20, 0.30, 0.50),
     right=0.99, left=0.1,
     top=0.93, bottom=0.08,
     wspace=0.15, hspace=0.40
)

col1_ylabel = -0.06
col2_ylabel = -0.22
#% =====================================================Make panel A===================================================== %%#
nlims = (-0.1, 1.1)
vlims = (-55, 5)

#Do the plotting
gsAL = gs[1, 1].subgridspec(ncols=1, nrows=2)
axA1 = fig1.add_subplot(gsAL[1])
axA1.plot(tSpike, vSpike, c=v_color, lw=3.0)
ylabel("Vt")
axA1.yaxis.set_label_coords(col1_ylabel, 0.5)
axA1.xaxis.set_visible(false) #Turn off the bottom axis

axA2 = fig1.add_subplot(gsAL[2])
axA2.plot(tSpike, nSpike, c=n_color, lw=3.0)
ylabel("Nt")
xlabel("Time (ms)")
axA2.yaxis.set_label_coords(col1_ylabel, 0.5)

axAR = fig1.add_subplot(gs[1, 2])
xlim(nlims)
ylim(vlims)
axAR.plot(nSpike, vSpike, c=:black, lw=3.0)
ylabel("Vt")
xlabel("Nt")
axAR.yaxis.set_label_coords(col2_ylabel, 0.5)

#% ===============================================Make panel B=============================================== %%#
gsBL = gs[2, 1].subgridspec(ncols=1, nrows=2)
axB1 = fig1.add_subplot(gsBL[1])
axB1.plot(tBurst, vBurst, c=v_color, lw=3.0)
axB1.xaxis.set_visible(false) #Turn off the bottom axis
ylabel("Vt (mV)")
axB1.yaxis.set_label_coords(col1_ylabel, 0.5)

axB2 = fig1.add_subplot(gsBL[2])
axB2.plot(tBurst, cBurst, c=c_color, lw=3.0)
xlabel("Time (s)")
ylabel("Ct")
axB2.yaxis.set_label_coords(col1_ylabel, 0.5)

axBR = fig1.add_subplot(gs[2, 2])
axBR.plot(vrng, crng, c=c_color, lw=3.0)
ylabel("Ct")
xlabel("Vt")
axBR.yaxis.set_label_coords(col2_ylabel, 0.5)

#% ===================================================Figure Panel C=================================================== %%#
gsCL = gs[3, 1].subgridspec(ncols=1, nrows=4)
axCL1 = fig1.add_subplot(gsCL[1])
axCL1.plot(tIBI, cIBI, c=c_color, lw=3.0)
ylabel("Ct")
axCL1.xaxis.set_visible(false) #Turn off the bottom axis
axCL1.yaxis.set_label_coords(col1_ylabel, 0.5)

axCL2 = fig1.add_subplot(gsCL[2])
axCL2.plot(tIBI, aIBI, c=a_color, lw=3.0)
ylabel("At")
axCL2.xaxis.set_visible(false) #Turn off the bottom axis
axCL2.yaxis.set_label_coords(col1_ylabel, 0.5)

axCL3 = fig1.add_subplot(gsCL[3])
axCL3.plot(tIBI, bIBI, c=b_color, lw=3.0)
ylabel("Bt")
axCL3.xaxis.set_visible(false) #Turn off the bottom axis
axCL3.yaxis.set_label_coords(col1_ylabel, 0.5)

axCL4 = fig1.add_subplot(gsCL[4])
axCL4.plot(tIBI, vIBI, c=v_color, lw=3.0)
ylabel("Vt")
xlabel("Time (s)")
axCL4.yaxis.set_label_coords(col1_ylabel, 0.5)

gsCR = gs[3, 2].subgridspec(ncols=1, nrows=5, height_ratios=[0.25, 0.125, 0.25, 0.125, 0.25])
axCR1 = fig1.add_subplot(gsCR[1])
axCR1.plot(aIBI, cIBI, c=a_color, lw=3.0)
ylabel("Ct")
xlabel("At")
axCR1.yaxis.set_label_coords(col2_ylabel, 0.5)

axCR2 = fig1.add_subplot(gsCR[3])
axCR2.plot(bIBI, aIBI, c=b_color, lw=3.0)
ylabel("At")
xlabel("Bt")
axCR2.yaxis.set_label_coords(col2_ylabel, 0.5)

axCR3 = fig1.add_subplot(gsCR[5])
axCR3.plot(vIBI, bIBI, c=v_color, lw=3.0)
ylabel("Bt")
xlabel("Vt")
axCR3.yaxis.set_label_coords(col2_ylabel, 0.5)
println(" Complete")

#%% Save the figure
loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2021 A Computational Model - Sci. Rep\Figures"
print("[$(now())]: Saving the figure 1...")
fig1.savefig("$(loc)/figure1_ModelVariables.png")
plt.close("all")
println(" Completed")