#%% This shows you how to run a Codimensional analysis
using Revise
using RetinalChaos
include("../figures/figure_setup.jl")

#%% Part 1: Running an ensemble problems
#Setp 1: Import all parameters and make the model
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict #Initial conditions
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 10.0
pars_dict[:ρe] = 0.0
pars_dict[:ρi] = 0.0
#pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict #Parameters
tspan = (0.0, 120e3) #Timespan
prob = ODEProblem(T_ODE, u0, tspan, p) #ODE problem

#Step 2: Determine the number of trajectories and the parameter to adjust
n_trajectories = 40
par_idx = p_find(:δ; list_p=t_pars); #Point to the index of the parameter
test_rng = LinRange(0.01, 0.05, n_trajectories); #Determine the range of the parameters (specified above)

#Step 3: Set up the ensemble problem
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #Set up the problem

#Step 4: Run the simulation
@time sim = solve(ensemble_prob, saveat=0.10, trajectories=n_trajectories, EnsembleThreads());

#[OPTIONAL]: Plot the solutions 
plt_a = plot(sim[1], vars=[1], c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
plt_b = plot(sim[1], vars=(1, 2), c=:jet, line_z=1, xlims=(-70.0, 10.0), clims=(test_rng[1], test_rng[end]))
for (sol_idx, sol_i) in enumerate(sim)
    println(test_rng[sol_idx])
    plot!(plt_a, sol_i, vars=[1], c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
    plot!(plt_b, sol_i, vars=(1, 2), c=:jet, xlims=(-70.0, 10.0), line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
end
plot(plt_a, plt_b, layout=2)

#%% Part 2: Running a equilibria analysis
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 10.0 #Set initial applied current to 0
pars_dict[:ρi] = 0.0 #remove GABA influence
pars_dict[:ρe] = 0.0 #remove ACh influence
pars_dict[:g_K] = 3.2
pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict

#Step 3: determine the timespan
tspan = (0.0, 300e3);

#Step 4: set up the equilibria problem
prob_eq = ODEProblem(T_ODE, u0, tspan, p)
sol = solve(prob_eq, progress=true, progress_steps=1)
plot(sol, vars=(1, 2), plotdensity=Int64(1e6))

#%% Conduct the equilibria analysis
@time eq_analysis = find_equilibria(prob_eq)
print(eq_analysis)
plot(eq_analysis)

@time eq_analysis = find_equilibria(prob_eq, equilibrium_resolution = 9)
plot(eq_analysis)


print(eq_analysis)
plot(eq_analysis)
#%% Step 5: Running a Codimensional-1 analysis
codim1 = (:I_app)
print("Codimensional analysis over parameter $codim1")
c1_lims = (-10.0, 200.0)
print("On parameter range: $c1_lims beginning:")
@time c1_map = codim_map(prob_eq, codim1, c1_lims, equilibrium_resolution=10)
println("Complete")

#%% Plot the bifurcation analysis
eq_plot = plot(c1_map, xlabel="Injected Current", ylabel="Membrane Voltage")
bif_val, bif_eq = find_bifurcation(c1_map)
saddle_vs = map(x -> x.saddle[1][1], bif_eq)
plot!(eq_plot, bif_val, saddle_vs, marker=:square, c=:blue, seriestype=:scatter, label="Saddle node bif")


#%% Step 6: Run out all of the solutions contained in the codimensional analysis
prob_eq = ODEProblem(T_ODE, u0, tspan, p) #First we need to reset the function
test_rng = map(x -> x[1], c1_map.points) #Determine the range of the parameters (specified above)
par_idx = codim1 |> p_find #Point to the index of the parameter
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob_eq, prob_func=prob_func); #Set up the problem
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, trajectories=length(test_rng), EnsembleThreads(), progress=true, progress_steps=1);

#%% Plot Codim Solutions
for (sol_idx, sol) in enumerate(sim)
    dt = 100.0
    vt = map(t -> sol(t)[1], collect(5e3:dt:tspan[2]))
    zt = repeat([test_rng[sol_idx]], length(vt))
    plot!(eq_plot, zt, vt, legend=false, marker=:circle, c=:blue)
end
eq_plot

#%% Step 6 Codim 2 analysis
codim2 = (:g_K, :g_Ca)
c1_lims = (1.0, 5.0);
c2_lims = (5.0, 15.0);
print("Codimensional analysis time to complete:")
@time c2_map = codim_map(prob_eq, codim2, c1_lims, c2_lims);
#%% Find the 2D bifurcation
find_bifurcation(c2_map)
#%%
plot(c2_map, view=:xyz, xlabel="I_app", ylabel="g_Ca", zlabel="Equilibrium Volt", legend=true)

#%% Codim 2 analysis
codim2 = (:g_K, :I_app)
c1_lims = (0.0, 20.0);
c2_lims = (-50.0, 1.0);
print("Codimensional analysis time to complete:")
@time c2_map = codim_map(prob_eq, codim2, c1_lims=c1_lims, c2_lims=c2_lims);
plot(c2_map, view=:yx, xlabel="I_app", ylabel="g_K", legend=true)

#%% Codim 2 analysis between gCa and gK
codim2 = (:g_Ca, :g_K)
c1_lims = (5.0, 30.0);
c2_lims = (1.0, 20.0);
print("Codimensional analysis time to complete:")
@time c2_map = codim_map(prob_eq, codim2, c1_lims=c1_lims, c2_lims=c2_lims);
#%%
plt = plot(c2_map, view=:xy, xlabel="g_Ca", ylabel="g_K")#, #xlabel = "I_app", ylabel = "g_Ca", legend = true)
#%%
savefig(plt, "figures\\supp_dyn.png")
#%%
c3_map
#%% Testing noise plots
p = read_JSON(params_file);
p[:σ] = 0.01
p[:τw] = 1000
u0 = read_JSON(conds_file);
tspan = (0.0, 300e3)
prob = SDEProblem(T_sde, u0 |> extract_dict, tspan, p |> extract_dict);
#sol = solve(prob)
#plot(sol, vars = [:v])
#%% Iterate through the 
print("Time it took to simulate $(tspan[2]/1000)s:")
test_rng = range(1.0, 30.0, length=50) #this ranges from halving the parameter to doubling it
par_idx = findall(x -> x == :g_Leak, Symbol.(T_sde.ps))
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng)
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func);
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, trajectories=length(test_rng), EnsembleThreads());
trace_plot = plot(legend=false)
for (sol_idx, sol_i) in enumerate(sim)
    plot!(trace_plot, sol_i, vars=[:v], line_z=test_rng[sol_idx], c=:jet, zlim=(1, length(sim)), layout=grid(2, 1), colorbar=true)
end
results = zeros(3, length(test_rng))
for (sol_idx, sol) in enumerate(sim)
    dt = 0.1 #set the time differential according to supp figure 1
    t_rng = collect(sol.t[1]:dt:sol.t[end]) #set the time range
    v_t = map(t -> sol(t)[1], t_rng) #extract according to the interval
    ts_analysis = timescale_analysis(v_t, dt=dt)
    for i = 1:length(ts_analysis)
        results[i, sol_idx] = sum(ts_analysis[i]) / length(ts_analysis[i])
    end
end
p1 = plot(test_rng, results[1, :])
p2 = plot(test_rng, results[2, :])
p3 = plot(test_rng, results[3, :])
sfig1 = plot(trace_plot, p1, p2, p3, layout=grid(4, 1, heights=[0.70, 0.10, 0.10, 0.10]))

#%% Testing synchrony of bursts
p = read_JSON(params_file);
p[:σ] = 0.25
p[:g_ACh] = 0.0
u0 = read_JSON(conds_file);
tspan = (0.0, 120e3)
prob = SDEProblem(T_sde, u0 |> extract_dict, tspan, p |> extract_dict);
print("Time it took to simulate $(tspan[2]/1000)s:")
sol = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, progress=true);
plot(sol, vars=[:v, :e], layout=grid(2, 1))