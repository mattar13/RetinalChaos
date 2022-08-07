using RetinalChaos
import RetinalChaos: calculate_threshold, get
import RetinalChaos: extract_equilibria, find_equilibria
using LaTeXStrings, Statistics
include("figure_setup.jl");

println("Running the plotting script for figure 2")

#%% Open data
print("[$(now())]: Setting up modelling data... ")
#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
conds_dict[:v] = -25.0
u0 = conds_dict |> extract_dict
#Step 2: Import the parameters
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 15.0
pars_dict[:ρi] = 0.0
pars_dict[:ρe] = 0.0
p = pars_dict |> extract_dict
#Step 3: determine the timespan
tspan = (0.0, 300e3)
#Step 4: set up the problem
prob = SDEProblem(T_SDE, noise, u0, tspan, p)
#Step 5: Solve the problem
sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1);
println(" Completed")


#% Run a dynamical analysis to get the equilibrium
print("[$(now())]: Setting up equilibrium data... ")
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
pars_dict_eq = read_JSON("params\\params.json")
pars_dict_eq[:I_app] = 0.0 #Set initial applied current to 0
pars_dict_eq[:ρi] = 0.0 #remove GABA influence
pars_dict_eq[:ρe] = 0.0 #remove ACh influence
pars_dict_eq[:g_TREK] = 0.0 #Remove the sAHP
p_eq = pars_dict_eq |> extract_dict
prob_eq = ODEProblem(T_ODE, u0, (0.0, 100.0), p_eq)
# Conduct the codim analysis
codim1 = (:I_app)
c1_lims = (-70.0, 30.0)
c1_map = codim_map(prob_eq, codim1, c1_lims, equilibrium_resolution=10)
println(" Completed")
#% 
print("[$(now())]: Extracting data... ")
res = extract_equilibria(c1_map) #Pass back all of the equilibria
points = res[1]
saddle_p = res[2]
stable_p = res[3]
unstable_p = res[4]
unstable_focus_p = res[5]
stable_focus_p = res[6]

last_saddle_idx = findlast(isnan.(saddle_p) .== 0)
saddle_bifurcation = points[last_saddle_idx]
saddle_eq = saddle_p[last_saddle_idx]

#% Extract plotting data
dt = 1.0
t = (sol.t[1]:dt:120e3)
vt = sol(t, idxs=1)
bt = sol(t, idxs=5)
wt = sol(t, idxs=8)
t = t / 1000

gTREK = -pars_dict[:g_TREK]
EK = pars_dict[:E_K]
ITREK = gTREK .* bt[2:end] .* (vt[1:end-1] .- EK) #We want to measure the current from the value before
ITREK_avg = sum(ITREK)/length(ITREK)
ITREK_SEM = std(ITREK)/sqrt(length(ITREK))
println(" Completed")

#%% Let s set up the figures
print("[$(now())]: Plotting... ")
width_inches = 7.5
height_inches = 5.0
fig2 = plt.figure("Biophysical Noise", figsize=(width_inches, height_inches))

gs = fig2.add_gridspec(3, 2,
     width_ratios=(0.55, 0.45),
     height_ratios=(0.30, 0.30, 0.40),
     right=0.99, left=0.14,
     top=0.91, bottom=0.1,
     wspace=0.4, hspace=0.40
)
col1_ylabel = -0.15
col2_ylabel = -0.17

#% =====================================================Make panel A===================================================== %%#
axA = fig2.add_subplot(gs[1, 1])
ylim(-100.0, 0.0)
axA.plot(t, vt, c=v_color, lw=lw_standard)
ylabel("Membrane \n Voltage (mV)")
axA.xaxis.set_visible(false) #Turn off the bottom axis
axA.yaxis.set_label_coords(col1_ylabel, 0.5)
axA.yaxis.set_major_locator(MultipleLocator(50.0))
axA.yaxis.set_minor_locator(MultipleLocator(25.0))
axA.spines["bottom"].set_visible(false)
# Plot the second column
axA2 = fig2.add_subplot(gs[1, 2])
xlim(c1_lims)
ylim(-100.0, 0.0)
axA2.plot(points, saddle_p, c=:blue, lw=lw_standard)
axA2.plot(points, stable_p, c=:green, lw=lw_standard)
axA2.plot(points, stable_focus_p, c=:green, ls="--", lw=lw_standard)
axA2.plot(saddle_bifurcation, saddle_eq, marker="s", markersize=5.0, c=:cyan)
I_sn_text = L"I_{sn}"
axA2.legend(["Saddle Eq.", "Stable Eq.", "Stable Oscillation", "Bifurcation ($(I_sn_text))"],
     ncol=2, columnspacing=0.70,
     bbox_to_anchor=(0.04, 1.3), fontsize=9.0, markerscale=0.5, handletextpad=0.3
)
ylabel("Equilibria \n Voltage (mV)")
xlabel("Injected current")
#Plot the points where the saddle node dissappears
axA2.xaxis.set_visible(false) #Turn off the bottom axis
axA2.yaxis.set_label_coords(col2_ylabel, 0.5)
axA2.yaxis.set_major_locator(MultipleLocator(50.0))
axA2.yaxis.set_minor_locator(MultipleLocator(25.0))
axA2.spines["bottom"].set_visible(false)
#% ===============================================Make panel B=============================================== %%#
axB = fig2.add_subplot(gs[2, 1])
ylim(-60, 0.0)
axB.plot(t[1:end-1], ITREK, c=b_color, lw=lw_standard)
ylabel("TREK \n current (pA)")
axB.xaxis.set_visible(false) #Turn off the bottom axis
axB.yaxis.set_label_coords(col1_ylabel, 0.5)
axB.yaxis.set_major_locator(MultipleLocator(30.0))
axB.yaxis.set_minor_locator(MultipleLocator(15.0))
axB.spines["bottom"].set_visible(false)

axB2 = fig2.add_subplot(gs[2, 2])
xlim(c1_lims)
ylim(-100.0, 0.0)
axB2.plot(ITREK, vt[1:end-1], linewidth=1.0, c=b_color, marker="o", markersize=4.0, markerfacecolor=:black)
axB2.legend(["TREK current"],
     bbox_to_anchor=(0.04, 1.2), fontsize=9.0, handletextpad=0.3
)

axB2.plot(points, saddle_p, c=:blue, lw=lw_standard)
axB2.plot(points, stable_p, c=:green, lw=lw_standard)
axB2.plot(points, stable_focus_p, c=:green, ls="--", lw=lw_standard)
axB2.plot(saddle_bifurcation, saddle_eq, marker="s", markersize=5.0, c=:cyan)
#Plot the points where the saddle node dissappears
ylabel("Membrane \n Voltage (mV)")
xlabel("TREK \n Current (pA)")
axB2.xaxis.set_visible(false) #Turn off the bottom axis
axB2.yaxis.set_label_coords(col2_ylabel, 0.5)
axB2.yaxis.set_major_locator(MultipleLocator(50.0))
axB2.yaxis.set_minor_locator(MultipleLocator(25.0))
axB2.spines["bottom"].set_visible(false)
#% ===================================================Figure Panel C=================================================== %%#
axC = fig2.add_subplot(gs[3, 1])
ylim(-8.0, 8.0)
axC.plot(t, wt, c=:black, lw=lw_standard)
ylabel("Noisy Current \n (pA)")
xlabel("Time (s)")
axC.yaxis.set_label_coords(col1_ylabel, 0.5)
axC.yaxis.set_major_locator(MultipleLocator(8.0))
axC.yaxis.set_minor_locator(MultipleLocator(4.0))

axC.xaxis.set_major_locator(MultipleLocator(20.0))
axC.xaxis.set_minor_locator(MultipleLocator(10.0))

axC2 = fig2.add_subplot(gs[3, 2])
xlim(c1_lims)
ylim(0.0, 1.0)
hfit = fit(Histogram, wt, LinRange(c1_lims[1], c1_lims[2], 200))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
axC2.plot(edges, weights, c=:black, lw=lw_standard)
axC2.fill_between(edges[edges.<=saddle_bifurcation], weights[edges.<=saddle_bifurcation], color=:gray)
axC2.fill_between(edges[edges.>saddle_bifurcation], weights[edges.>saddle_bifurcation], color=:green)
axC2.legend(["Noise Histogram", "Non-spiking noise", "Spiking Noise"],
     bbox_to_anchor=(0.04, 1.0), fontsize=9.0, handletextpad=0.3
)

ylabel("Probability of \n Noisy Current")
xlabel("Noise Current (pA)")
axC2.yaxis.set_label_coords(col2_ylabel, 0.5)
axC2.yaxis.set_major_locator(MultipleLocator(0.5))
axC2.yaxis.set_minor_locator(MultipleLocator(0.25))
axC2.xaxis.set_major_locator(MultipleLocator(20.0))
axC2.xaxis.set_minor_locator(MultipleLocator(10.0))
println(" Completed")

axC2.annotate("A", (0.01, 0.90), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axC2.annotate("B", (0.53, 0.90), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axC2.annotate("C", (0.01, 0.65), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axC2.annotate("D", (0.53, 0.65), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")
axC2.annotate("E", (0.01, 0.35), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

axC2.annotate("F", (0.53, 0.35), xycoords="figure fraction", annotation_clip=false, fontsize=20.0, fontweight="bold")

#%% Save the Plot 
loc = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Journal Submissions\2022 A Computational Model - Sci. Rep\Submission 1\Figures"
print("[$(now())]: Saving the figure 2...")
fig2.savefig("$(loc)/figure2_BiophysicalProperties.jpg")
plt.close("all")
println(" Completed")