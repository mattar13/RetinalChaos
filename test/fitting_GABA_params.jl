using Revise
using RetinalChaos #Load the package
using LinearAlgebra
using StaticArrays
using Plots
#%%Lets figure out how to extract the waves
import RetinalChaos: GABA_conds, GABA_pars, T_ODE, GABA_ODE, GABA_SDE

#%% Lets run a baseline experiment to compare traces
conds_dict_GABA = read_JSON("params/GABA_conds.json")
u0 = extract_dict(conds_dict_GABA, GABA_conds)

#Step 2: Import the parameters
pars_dict = read_JSON("params/GABA_params.json")
pars_dict[:g_GABA] = 15.0
p_GABA = extract_dict(pars_dict, GABA_pars)

pars_dict[:g_GABA] = 0.0 #Basically removes the influence of the GABA receptors
p = extract_dict(pars_dict, GABA_pars)
#Step 3: determine the timespan
tspan = (0.0, 300e3);

#Step 4 Set up the two problems
prob = ODEProblem(GABA_ODE, u0, tspan, p); #1 is the Plain ODE problem
probGABA = ODEProblem(GABA_ODE, u0, tspan, p_GABA); #1 is the Plain ODE problem

@time sol = solve(prob, progress=true);
@time solGABA = solve(probGABA, progress=true);

# Plotting results
plt_a = plot(sol, vars=[1, 6, 7], layout=(3, 1))
plot!(plt_a, solGABA, vars=[1, 6, 7], layout=(3, 1))


#%% Lets run a range of different g_GABA params
n = 10
probGABA = ODEProblem(GABA_ODE, u0, tspan, p_GABA); #1 is the Plain ODE problem
par_idx = p_find(:g_GABA; list_p=GABA_pars) #Point to the index of the parameter
test_rng = LinRange(5.0, 20.0, n) #Determine the range of the parameters (specified above)
prob_func(prob, i, repeat) = ensemble_func(probGABA, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(probGABA, prob_func=prob_func); #Set up the problem
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, trajectories=n, EnsembleThreads());

plt_a = plot(sim[1], vars=[1], c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
plt_b = plot(sim[1], vars=(1, 2), c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
for (sol_idx, sol_i) in enumerate(sim)
     println(test_rng[sol_idx])
     plot!(plt_a, sol_i, vars=[1], c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
     plot!(plt_b, sol_i, vars=(1, 2), c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
end
plot(plt_a, plt_b, layout=2)

#%% Lets run this as a Stochastic problem
conds_dict_GABA = read_JSON("params/GABA_conds.json")
u0 = extract_dict(conds_dict_GABA, GABA_conds)

pars_dict = read_JSON("params/GABA_params.json")
pars_dict[:g_GABA] = 13
p_GABA = extract_dict(pars_dict, GABA_pars)
pars_dict[:g_GABA] = 0.0 #Basically removes the influence of the GABA receptors
p = extract_dict(pars_dict, GABA_pars)
tspan = (0.0, 300e3);

probSDE = SDEProblem(GABA_SDE, noise, u0, tspan, p); #1 is the Plain ODE problem
probSDE_GABA = SDEProblem(GABA_SDE, noise, u0, tspan, p_GABA); #1 is the Plain ODE problem

@time solSDE = solve(probSDE, SOSRI());
@time solSDE_GABA = solve(probSDE_GABA, SOSRI());

plt_a = plot(solSDE, vars=[1, 6, 7, 8], layout=(4, 1))
plot!(plt_a, solSDE_GABA, vars=[1, 6, 7, 8], layout=(4, 1))


#%% We need to adjust the diffusion stencils
import RetinalChaos.GABA_PDE
nx = ny = 50;
p_dict = read_JSON(params_file) #set up parameters
p_dict[:g_GABA] = 1.2
u_dict = read_JSON(conds_file) #Set up the initial conditions
p_net = extract_dict(p_dict, GABA_pars)
u0_net = extract_dict(u_dict, GABA_conds, (nx, ny))
tspan = (0.0, 120e3)
probPDE = SDEProblem(GABA_PDE, noise, u0_net, tspan, p_net)
@time NetSol = solve(probPDE, SOSRI(),
     abstol=2e-2, reltol=0.2, maxiters=1e7,
     save_idxs=[1:(nx*ny)...],
     progress=true, progress_steps=1
)

#%%
dFrame = 100.0
anim = @animate for t in NetSol.t[1]:dFrame:NetSol.t[end]
     println(t)
     frame = reshape(NetSol(t), (nx, ny))
     plot(frame, st=:heatmap, ratio=:equal, c=:jet, clims=(-70.0, 10.0))
end

gif(anim, fps=10.0)

plot(Array(NetSol)[2, :])