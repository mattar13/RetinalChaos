#%% This shows you how to run a Codimensional analysis
#using Revise
using RetinalChaos
include("../figures/figure_setup.jl")
import BifurcationKit.@set
#reload_parameters()
#Step 1: Import:
import RetinalChaos.SDEModel #the SDEModel
import RetinalChaos.ODEModel #the ODEModel
import RetinalChaos.u0 #the initial condutions 
import RetinalChaos.parameters #the parameters

#Step 2: determine the timespan
tmin = 0.0
tmax = 500.0

#%%=Run equilibrium for chloride==========================================================================================================#
#Step 3: set up the problem making sure enable jacobian is true
reload_parameters()
parameters[I_app] = 25.0
parameters[g_ACh] = 0.0
#parameters[g_GABA] = 0.0
parameters[g_TREK] = 0.0 #Remove the sAHP
u0[i] = parameters[ρi]
u0[e] = parameters[ρe] = 0.0
@time ODEProb = ODEProblem(ODEModel, u0, (tmin, tmax), parameters, jac=true)
sol = solve(ODEProb)
initial_plot = plot(sol, idxs=[v, e, i], layout = 3)

# Part A: Determining jacobians and gradients
#Step 1: Make the function and Jacobian
ODEFunc = ODEProb.f
F = (u, p) -> ODEFunc(u, p, 0.0)
J = (u, p) -> ODEFunc.jac(u, p, 0.0)

id_par = indexof(E_Cl)
par_tm = ODEProb.p

prob = BK.BifurcationProblem(
    F, ODEProb.u0, par_tm, (@lens _[id_par]);
    J=J,
    recordFromSolution=(x, p) -> (V=x[2], N=x[3])
)

# continuation options (Theses need to be tuned)
opts_br = BK.ContinuationPar(
    pMin=-80.0, pMax=-40.0, # parameters to have a smooth result
    ds=0.04, dsmax=1.0,# this is to detect bifurcation points precisely with bisection
    detectBifurcation=3, # Optional: bisection options for locating bifurcations
    nInversion=8, maxBisectionSteps=50,
    maxSteps=500,
    nev=3
)

continuation_method = BK.PALC(tangent=BK.Bordered())

br = BK.continuation(prob, continuation_method, opts_br;
    normC=norminf,
    verbosity=3,
    bothside=true, plot=true
)

# This section determines the periodic orbits 
optn_po = BK.NewtonPar(verbose=true, tol=1e-8, maxIter=100)
opts_po_cont = BK.ContinuationPar(
    dsmin=0.04, ds=0.1, dsmax=2.0,
    pMin=0.0, pMax=50.0,
    maxSteps=110,
    newtonOptions=(BK.@set optn_po.tol = 1e-7),
    nev=3, plotEveryStep=2, detectBifurcation=0,
    saveSolEveryStep=1
)

#This is for printing
args_po = (recordFromSolution=(x, p) -> begin
        xtt = BK.getPeriodicOrbit(p.prob, x, @set par_tm[id_par] = p.p)
        return (max=maximum(xtt[1, :]),
            min=minimum(xtt[1, :]),
            period=BK.getPeriod(p.prob, x, @set par_tm[id_par] = p.p))
    end,
    plotSolution=(x, p; k...) -> begin
        xtt = BK.getPeriodicOrbit(p.prob, x, @set par_tm[id_par] = p.p)
        plot!(xtt.t, xtt[2, :]; label="V", k...)
        plot!(xtt.t, xtt[3, :]; label="N", k...)
        plot!(br; subplot=1, putspecialptlegend=false)
    end,
    normC=norminf
)

Mt = 15 # number of time sections
@time br_pocoll = BK.continuation(
    br, 4, opts_po_cont, # we want to branch form the 4th bif. point
    BK.PeriodicOrbitOCollProblem(Mt, 5), # we want to use the Collocation method to locate PO, with polynomial degree 5
    plot=true, verbosity=3;
    args_po... # regular continuation options
)
scene = plot(br, br_pocoll, markersize=3)
plot!(scene, br_pocoll.param, br_pocoll.min, label="")

#periodic_orbits = plot()
# fetch the saved solutions
for sol in br_pocoll.sol
    # periodic orbit
    po = sol.x
    # get the mesh and trajectory
    traj = BK.getPeriodicOrbit(br_pocoll.prob, po, @set par_tm[id_par] = sol.p)
    plot!(scene, traj[1, :], traj[2, :], xlabel="E", ylabel="x", label="")
end
scene

#%% Codimensional analysis


#Step 4: set the parameters
reload_parameters()
parameters[E_Cl]
prob = SDEProblem(SDEModel, u0, (tmin, tmax), parameters) #ODE problem
#prob = SDEProblem(T_SDE, noise, u0, tspan, p) #ODE problem

#Step 2: Determine the number of trajectories and the parameter to adjust
reload_parameters()
n_trajectories = 5
par = :E_Cl #Adjust the noise
pmin = -65.0 #Determine the minimum parameter
pmax = 55.0 #Determine the maximum range of parameters
test_rng = LinRange(pmin, pmax, n_trajectories); #Determine the range of the parameters (specified above)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func); #Set up the problem

#Step 4: Run the simulation #Noise uses SOSRI(), 
@time sim = solve(ensemble_prob, SOSRI(), saveat=1.0, trajectories=n_trajectories, EnsembleThreads());

#%% Plot the solutions 
plt_a = plot(sim[1], idxs=[v], c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
plt_b = plot(sim[1], idxs=(v, n), c=:jet, line_z=1, ylims=(0.0, 1.0), xlims=(-70.0, 10.0), clims=(test_rng[1], test_rng[end]))
for (sol_idx, sol_i) in enumerate(sim)
    println(test_rng[sol_idx])
    plot!(plt_a, sol_i, idxs=[1], c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
    plot!(plt_b, sol_i, idxs=(1, 2), c=:jet, ylims=(0.0, 1.0), xlims=(-70.0, 10.0), line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
end
plot(plt_a, plt_b, layout=2)

#%% Part 2: Using bifurcation kits
using Bifuc


parameters[I_app] = 10.0 #Set initial applied current to 0
parameters[ρi] = 0.0 #remove GABA influence
parameters[ρe] = 0.0 #remove ACh influence
parameters[g_TREK] = 0.0 #Remove the sAHP

#Step 3: determine the timespan
tmin = 0.0
tmax = 300e3

#Step 4: set up the equilibria problem
prob_eq = ODEProblem(ODEModel, u0, (tmin, tmax), parameters)
sol = solve(prob_eq, progress=true, progress_steps=1)
plot(sol, vars=(v, n), plotdensity=Int64(1e6))

#%% Step 5: Running a Codimensional-1 analysis
codim1 = (:δ)
print("Codimensional analysis over parameter $codim1")
c1_lims = (0.001, 0.10)
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