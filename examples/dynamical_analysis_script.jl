#%% This shows all the different things you can do with the dynamical analysis of the model
using RetinalChaos
#Set up the plotting
font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)

#Set up the file root and default parameters
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

#save everything in the figures folder
save_figs = "figures\\"
if isdir(save_figs) == false
    #The directory does not exist, we have to make it 
    mkdir(save_figs)
end

#%% Running a plain simulations
p = read_JSON(params_file);
p[:σ] = 0.01
p[:τw] = 1000
u0 = read_JSON(conds_file);
tspan = (0.0, 600e3) #10 mins long
prob = SDEProblem(T_sde, u0|>extract_dict, tspan, p|>extract_dict);
sol = solve(prob; progress = true)

plot(sol, vars = [:v])
#%% Supplemental figure, Bifurcation analysis of voltage injections
p = read_JSON(params_file) 
u0 = read_JSON(conds_file)
p[:I_app] = 0.0 #Set initial applied current to 0
p[:g_ACh] = 0.0 #Remove g_ACh influence
p[:g_TREK] = 0.0 #Remove g_TREK influence
tspan = (0.0, 30e3);
prob_eq = ODEProblem(T_ode, u0|>extract_dict, tspan, p|>extract_dict)
#We can superimpose a SDE problem overtop
prob_sde = SDEProblem(T_sde, u0|>extract_dict, tspan, p|>extract_dict);
#%% Conducting an Equilibrium analysis
sol_eq = find_equilibria(prob_eq)
#%%
sol_sde = solve(prob_sde)
#%% Codim 1 analysis
codim1 = (:I_app)
c1_lims = (-50.0, 50.0)
print("Codimensional analysis time to complete:")
@time c1_map = codim_map(prob_eq, codim1, c1_lims = c1_lims, eq_res = 10)
eq_plot = plot(c1_map, xlabel = "Injected Current", ylabel = "Membrane Voltage")
#Run a ensemble function (after resetting the function)
prob_eq = ODEProblem(T_ode, u0|>extract_dict, tspan, p|>extract_dict)
test_rng = range(c1_lims[1], c1_lims[2], length = 25) #this ranges from halving the parameter to doubling it
par_idx = findall(x -> x==codim1, Symbol.(T_sde.ps))
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng)
ensemble_prob = EnsembleProblem(prob_eq, prob_func = prob_func);
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, trajectories = length(test_rng), EnsembleThreads());
for (sol_idx, sol) in enumerate(sim)
    dt = 100.0
    vt = map(t -> sol(t)[1], collect(5e3:dt:tspan[2]))
    zt = repeat([test_rng[sol_idx]], length(vt))
    plot!(eq_plot, zt, vt, legend = false, marker = :circle,  c = :blue)
end
eq_plot

#%% Codim 2 analysis
codim2 = (:g_Ca, :I_app)
c1_lims = (0.0, 20.0); c2_lims = (-50.0, 1.0) 
print("Codimensional analysis time to complete:")
@time c2_map = codim_map(prob_eq, codim2, c1_lims = c1_lims, c2_lims = c2_lims);
plot(c2_map, view = :yx, xlabel = "I_app", ylabel = "g_Ca", legend = true)

#%% Codim 2 analysis
codim2 = (:g_K, :I_app)
c1_lims = (0.0, 20.0); c2_lims = (-50.0, 1.0) 
print("Codimensional analysis time to complete:")
@time c2_map = codim_map(prob_eq, codim2, c1_lims = c1_lims, c2_lims = c2_lims);
plot(c2_map, view = :yx, xlabel = "I_app", ylabel = "g_K", legend = true)

#%% Codim 2 analysis between gCa and gK
codim2 = (:g_Ca, :g_K)
c1_lims = (5.0, 30.0); c2_lims = (1.0, 20.0);
print("Codimensional analysis time to complete:")
@time c2_map = codim_map(prob_eq, codim2, c1_lims = c1_lims, c2_lims = c2_lims);
#%%
plt = plot(c2_map, view = :xy, xlabel = "g_Ca", ylabel = "g_K")#, #xlabel = "I_app", ylabel = "g_Ca", legend = true)
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
prob = SDEProblem(T_sde, u0|>extract_dict, tspan, p|>extract_dict);
#sol = solve(prob)
#plot(sol, vars = [:v])
#%% Iterate through the 
print("Time it took to simulate $(tspan[2]/1000)s:")
test_rng = range(1.0, 30.0, length = 50) #this ranges from halving the parameter to doubling it
par_idx = findall(x -> x==:g_Leak, Symbol.(T_sde.ps))
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng)
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func);
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, SOSRI(), abstol = 2e-2, reltol = 2e-2, maxiters = 1e7, trajectories = length(test_rng), EnsembleThreads());
trace_plot = plot(legend = false)
for (sol_idx, sol_i) in enumerate(sim)
    plot!(trace_plot, sol_i, vars = [:v], line_z = test_rng[sol_idx], c = :jet, zlim = (1, length(sim)), layout = grid(2,1), colorbar = true)    
end
results = zeros(3, length(test_rng))
for (sol_idx, sol) in enumerate(sim)
    dt = 0.1 #set the time differential according to supp figure 1
    t_rng = collect(sol.t[1]:dt:sol.t[end]) #set the time range
    v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
    ts_analysis = timescale_analysis(v_t, dt = dt)
    for i = 1:length(ts_analysis)
        results[i, sol_idx] = sum(ts_analysis[i])/length(ts_analysis[i])
    end
end
p1 = plot(test_rng, results[1, :])
p2 = plot(test_rng, results[2, :])
p3 = plot(test_rng, results[3, :])
sfig1 = plot(trace_plot, p1, p2, p3, layout = grid(4,1, heights = [0.70, 0.10, 0.10, 0.10]))

#%% Testing synchrony of bursts
p = read_JSON(params_file);
p[:σ] = 0.25
p[:g_ACh] = 0.0
u0 = read_JSON(conds_file);
tspan = (0.0, 120e3)
prob = SDEProblem(T_sde, u0|>extract_dict, tspan, p|>extract_dict);
print("Time it took to simulate $(tspan[2]/1000)s:")
sol = solve(prob, SOSRI(), abstol = 2e-2, reltol = 2e-2, maxiters = 1e7, progress = true);
plot(sol, vars = [:v, :e], layout = grid(2,1))