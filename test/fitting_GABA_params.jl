using Revise
using RetinalChaos #Load the package
using LinearAlgebra
using Plots
import RetinalChaos: GABA_conds, GABA_pars, T_ODE, GABA_ODE, GABA_SDE

#%% Lets run a baseline experiment to compare traces
conds_dict_GABA = read_JSON("params/GABA_conds.json")
u0 = extract_dict(conds_dict_GABA, GABA_conds)

#Step 2: Import the parameters
pars_dict = read_JSON("params/GABA_params.json")
pars_dict[:g_GABA] = 10.5
p_GABA = extract_dict(pars_dict, GABA_pars)

pars_dict[:g_GABA] = 0.0 #Basically removes the influence of the GABA receptors
p = extract_dict(pars_dict, GABA_pars)
#Step 3: determine the timespan
tspan = (0.0, 300e3);

#Step 4 Set up the two problems
prob = ODEProblem(GABA_ODE, u0, tspan, p); #1 is the Plain ODE problem
probGABA = ODEProblem(GABA_ODE, u0, tspan, p_GABA); #2 is the GABA ODE problem
probSDE = SDEProblem(GABA_SDE, noise, u0, tspan, p); #3 is the Plain SDE problem
probSDE_GABA = SDEProblem(GABA_SDE, noise, u0, tspan, p_GABA); #4 is the GABA SDE problem

@time sol = solve(prob, progress=true);
@time solGABA = solve(probGABA, progress=true);
@time solSDE = solve(probSDE, SOSRI());
@time solSDE_GABA = solve(probSDE_GABA, SOSRI());

# Plotting results
plt_a = plot(sol, vars=[1, 8], layout=(2, 1))
plot!(plt_a, solGABA, vars=[1, 8], layout=(2, 1))
plt_b = plot(solSDE, vars=[1, 8], layout=(2, 1))
plot!(plt_b, solSDE_GABA, vars=[1, 8], layout=(2, 1))

plot(plt_a, plt_b)
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

#%% Now lets run through several parameters using the wave model
nx = ny = 50
b = Binomial(1, 0.75) #This means 75% of cells are sensitive to GABA
null = Array{Float64}(rand(b, nx, ny))
param = :g_GABA
for val in LinRange(8.0, 12.0, 10) #This is the range of values
     print("[$(now())]: Setting up binomal nullification... ")
     println(" [$(now())]: Completed")
     net = (dU, U, p, t) -> GABA_PDE_gNULL(dU, U, p, t, null)
     print("[$(now())]: Loading Parameters... ")
     conds_dict = read_JSON("params/GABA_conds.json")
     u0 = extract_dict(conds_dict, GABA_conds, dims=(nx, ny))
     pars_dict = read_JSON("params/GABA_params.json")
     #CHANGE ANY PARAMETER HERE
     pars_dict[param] = val
     p = extract_dict(pars_dict, GABA_pars)
     tspan_warmup = (0.0, 200e3)
     println(" [$(now())]: Completed")

     print("[$(now())]: Warming up model... ")
     prob = SDEProblem(net, noise, u0, tspan, p)
     @time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1)
     println(" [$(now())]: Completed")
     
     print("[$(now())]: Running model... ")
     prob = SDEProblem(net, noise, warmup[end], tspan, p)
     @time NetSol = solve(prob, SOSRI(), abstol=2e-2, reltol=0.2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...])
     println(" [$(now())]: Completed")

     # Step 7: Animate the solution
     print("[$(now())]: Animating solution... ")
     loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\GABA_experiments"
     animate_dt = 60.0
     anim = @animate for t = 1.0:animate_dt:NetSol.t[end]
          println("[$(now())]: Animating simulation...")
          frame_i = reshape(NetSol(t) |> Array, (nx, ny))
          heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-70.0, 0.0))
     end
     gif(anim, "$(loc)\\$(param)_$(round(val, digits = 1))animation.gif", fps=1000.0 / animate_dt)
     println(" [$(now())]: Completed")
end


#parameter spaces: 
#1) gGABA = [8.0->12.0]
#2) gCa = [10.0->20.0]
#3) gK = [7.5 -> 12.5]
#4) μ = [0.7-0.8]
