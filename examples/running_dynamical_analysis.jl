#%% This shows you how to run a Codimensional analysis
using Revise
using RetinalChaos

#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict

#Step 2: Import the parameters (and set some parameters to 0.0)
pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 0.0 #Set initial applied current to 0
pars_dict[:g_ACh] = 0.0 #Remove g_ACh influence
pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict

#Step 3: determine the timespan
tspan = (0.0, 120e3);

#Step 4: set up the equilibria problem
prob_eq = ODEProblem(T_ode, u0, tspan, p)
prob_sde = SDEProblem(T_sde, noise, u0, tspan, p);
eq_analysis = find_equilibria(prob_eq)
print(eq_analysis)

#%% Step 5: Running a Codimensional-1 analysis
codim1 = (:I_app)
print("Codimensional analysis over parameter $codim1")
c1_lims = (-60.0, 50.0)
print("On parameter range: $c1_lims beginning:")
@time c1_map = codim_map(prob_eq, codim1, c1_lims, equilibrium_resolution=10)
println("Complete")

#%% Step 6: Plot out all of the solutions contained in the codimensional analysis
prob_eq = ODEProblem(T_ode, u0, tspan, p) #First we need to reset the function
test_rng = map(x -> x[1], c1_map.points) #Determine the range of the parameters (specified above)
par_idx = p_find(codim1; list_p=RetinalChaos.tar_pars) #Point to the index of the parameter
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob_eq, prob_func=prob_func); #Set up the problem
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, trajectories=length(test_rng), EnsembleThreads());

#%% If plotting is loaded then plot the current vs voltage
eq_plot = plot(c1_map, xlabel="Injected Current", ylabel="Membrane Voltage")
bif_val, bif_eq = find_bifurcation(c1_map)
println(bif_val)
saddle_vs = map(x -> x.saddle[1][1], bif_eq)
plot!(eq_plot, bif_val, saddle_vs, marker=:square, c = :blue, seriestype=:scatter, label="Saddle node bif")

for (sol_idx, sol) in enumerate(sim)
    dt = 100.0
    vt = map(t -> sol(t)[1], collect(5e3:dt:tspan[2]))
    zt = repeat([test_rng[sol_idx]], length(vt))
    plot!(eq_plot, zt, vt, legend=false, marker=:circle, c=:blue)
end

eq_plot

#%% Step 6 Codim 2 analysis
codim2 = (:I_app, :g_Ca)
c1_lims = (-60.0, 50.0);
c2_lims = (0.0, 20.0);
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