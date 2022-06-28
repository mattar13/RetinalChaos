using Revise
using RetinalChaos
import RetinalChaos: calculate_threshold, get
import RetinalChaos: extract_equilibria, find_equilibria
include("figure_setup.jl");

#%% Open data
#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
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
prob = SDEProblem(T_SDE, noise, u0, tspan, p)
#Step 5: Solve the problem
@time sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true, progress_steps=1);

#%% Analyze the data
thresholds = calculate_threshold(sol)
spike_timestamps = get_timestamps(sol)
spike_durs, isi = extract_interval(spike_timestamps, max_duration=100, max_interval=100)
burst_tstamps, SPB = max_interval_algorithim(spike_timestamps)
burst_durs, ibi = extract_interval(burst_tstamps)

#%% Pull out example data (from a single burst)
dt = 1.0
t = sol.t[1]:dt:sol.t[end]
vt = sol(t, idxs=1)
bt = sol(t, idxs=5)
wt = sol(t, idxs=8)

#%% Run a dynamical analysis to get the equilibrium
conds_dict = read_JSON("params\\conds.json")
conds_dict[:v] = vm_noise
u0 = conds_dict |> extract_dict

pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = iapp_noise #Set initial applied current to 0
pars_dict[:ρi] = 0.0 #remove GABA influence
pars_dict[:ρe] = 0.0 #remove ACh influence
pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict
tspan = (0.0, 100.0)
prob_eq = ODEProblem(T_ODE, u0, tspan, p)
# Conduct the codim analysis
codim1 = (:I_app)
c1_lims = (45.0, 50.0)
@time c1_map = codim_map(prob_eq, codim1, c1_lims, equilibrium_resolution=10)

#%% Plot Codim Solutions
res = extract_equilibria(c1_map) #Pass back all of the equilibria
points = res[1]
saddle_p = res[2]
stable_p = res[3]
unstable_p = res[4]
unstable_focus_p = res[5]
stable_focus_p = res[6]
#bif_val, bif_eq = find_bifurcation(c1_map)
#saddle_vs = map(x -> x.saddle[1][1], bif_eq)
# Plot 
plot(points, saddle_p, c=:blue)
plot!(points, stable_p, c=:green)
plot!(points, unstable_p, c=:red)
plot!(points, stable_focus_p, c=:red, linestyle=:dash)
plot!(points, unstable_focus_p, c=:green, linestyle=:dash)


#%% Lets look deeper into some issues
eq_val = findall((isnan.(stable_p)) .== 0)[1]
i_app_eq = c1_map.points[eq_val][1]
u0_eq = c1_map.equilibria[eq_val].stable[1] #We will use this as out Initial condition

pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = i_app_eq #Set initial applied current to 0
pars_dict[:ρi] = 0.0 #remove GABA influence
pars_dict[:ρe] = 0.0 #remove ACh influence
pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict
tspan = (0.0, 100.0)
prob_eq = ODEProblem(T_ODE, u0_eq, tspan, p)
sol_eq = solve(prob_eq)
plot(sol_eq)

#%%
plt.clf()
#%% Let s set up the figures
width_inches = 16.0
height_inches = 10.0
fig2 = plt.figure("Biophysical Noise", figsize=(width_inches, height_inches))

gs = fig2.add_gridspec(3, 2,
     width_ratios=(0.70, 0.30),
     height_ratios=(0.30, 0.30, 0.40),
     right=0.99, left=0.1,
     top=0.93, bottom=0.08,
     wspace=0.15, hspace=0.40
)

axA = fig2.add_subplot(gs[1, 1])
axA.plot(t, vt, c=v_color, lw=3.0)
axB = fig2.add_subplot(gs[2, 1])
axB.plot(t, bt, c=b_color, lw=3.0)
axC = fig2.add_subplot(gs[3, 1])
axC.plot(t, wt, c=:black, lw=3.0)

#%% Plot the second column
axA2 = fig2.add_subplot(gs[1, 2])

axB2 = fig2.add_subplot(gs[2, 2])

axC2 = fig2.add_subplot(gs[3, 2])
hfit = fit(Histogram, wt, LinRange(-10, 10, 50))
weights = hfit.weights / maximum(hfit.weights)
edges = collect(hfit.edges[1])[1:length(hfit.weights)]
axC2.plot(weights, edges, c=:black, lw=3.0)
axC2.fill_between(weights, edges, color=:gray)
