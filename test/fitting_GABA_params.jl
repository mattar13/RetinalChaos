using Revise
using RetinalChaos #Load the package
using LinearAlgebra
using Plots
import RetinalChaos: GABA_conds, GABA_pars, T_ODE, GABA_ODE, GABA_SDE, ħ

#%% Lets replicate some GABA and ACh IV curves
v_rng = -120:20:100
pars_dict = read_JSON("params/GABA_params.json")
g_GABA = pars_dict[:g_GABA] = 0.5
g_ACh = pars_dict[:g_ACh] = 0.2
k_GABA = pars_dict[:k_GABA]
k_ACh = pars_dict[:k_ACh]
E_GABA = pars_dict[:E_GABA] = -60
E_ACh = pars_dict[:E_ACh]

I_GABA(v, i) = -g_GABA * ħ(i, k_GABA) * (v - E_GABA)
I_ACh(v, e) = -g_ACh * ħ(e, k_ACh) * (v - E_ACh)
loE = 0.0
hiE = 6.0
loI = 0.0
hiI = 5.0

plt_I = plot(label="GABA I-V")
plt_E = plot(label="ACh I-V")
plt_b = plot(label="Total I-V")
for nt in (0.0:0.1:5.0)
     plot!(plt_I, v_rng, v -> I_GABA(v, nt), c=:jet, line_z=nt, clims=(0.0, 5.0), label="")
     plot!(plt_E, v_rng, v -> I_ACh(v, nt), c=:jet, line_z=nt, clims=(0.0, 5.0), label="")
     plot!(plt_b, v_rng, v -> (I_ACh(v, nt) + I_GABA(v, nt)), c=:jet, line_z=nt, clims=(0.0, 5.0), label="")
end
plt_I
plt_E
plt_b
plot(plt_I, plt_E, plt_b, layout=(3, 1))


#%% Testing different Acetylcholine conditions
conds_dict_GABA = read_JSON("params/GABA_conds.json")
u0 = extract_dict(conds_dict_GABA, GABA_conds)
tspan = (0.0, 120e3)
#Step 2: Import the parameters
pars_dict = read_JSON("params/GABA_params.json")

pars_dict[:ρe] = 0.0
pars_dict[:ρi] = 5.0 #IN the presence of 5uM
p = extract_dict(pars_dict, GABA_pars)
prob_loE_hiI = SDEProblem(GABA_SDE, noise, u0, tspan, p); #1 is the Plain ODE problem

pars_dict[:ρi] = 0.0 #Low GABA
p = extract_dict(pars_dict, GABA_pars)
prob_loE_loI = SDEProblem(GABA_SDE, noise, u0, tspan, p); #1 is the Plain ODE problem

pars_dict[:ρe] = 6.0
pars_dict[:ρi] = 5.0 #IN the presence of 5uM
p = extract_dict(pars_dict, GABA_pars)
prob_hiE_hiI = SDEProblem(GABA_SDE, noise, u0, tspan, p); #1 is the Plain ODE problem

pars_dict[:ρi] = 0.0 #Low GABA
p = extract_dict(pars_dict, GABA_pars)
prob_hiE_loI = SDEProblem(GABA_SDE, noise, u0, tspan, p); #1 is the Plain ODE problem

@time sol_loE_hiI = solve(prob_loE_hiI, SOSRI(), progress=true);
@time sol_loE_loI = solve(prob_loE_loI, SOSRI(), progress=true);
@time sol_hiE_hiI = solve(prob_hiE_hiI, SOSRI(), progress=true);
@time sol_hiE_loI = solve(prob_hiE_loI, SOSRI(), progress=true);

# Plotting results
plt_a = plot(sol_loE_hiI, vars=[1], layout=(1, 1))
plot!(plt_a, sol_loE_loI, vars=[1], layout=(1, 1))
plt_b = plot(sol_hiE_hiI, vars=[1], layout=(1, 1))
plot!(plt_b, sol_hiE_loI, vars=[1], layout=(1, 1))

fig2_AChGABA = plot(plt_a, plt_b)

#%% Lets run a range of different g_GABA params
pars_dict = read_JSON("params/GABA_params.json")
#pars_dict[:ρe] = 0.0 #Eliminate ACh activity
p = extract_dict(pars_dict, GABA_pars)

n = 10
probGABA = SDEProblem(GABA_SDE, noise, u0, tspan, p); #1 is the Plain ODE problem
par_idx = p_find(:ρi; list_p=GABA_pars) #Point to the index of the parameter
test_rng = LinRange(0.0, 5.0, n) #Determine the range of the parameters (specified above)
prob_func(prob, i, repeat) = ensemble_func(probGABA, i, repeat, par_idx, test_rng) #Set up the problem function to indicate that the voltage will be altered
ensemble_prob = EnsembleProblem(probGABA, prob_func=prob_func); #Set up the problem
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, SOSRI(), trajectories=n, EnsembleThreads(), progress=true, progress_steps=1);

plt_a = plot(sim[1], layout=(2, 1), vars=[1, 7], c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
plt_b = plot(sim[1], vars=(1, 2), c=:jet, line_z=1, clims=(test_rng[1], test_rng[end]))
for (sol_idx, sol_i) in enumerate(sim)
     println(test_rng[sol_idx])
     plot!(plt_a, sol_i, vars=[1, 7], layout=2, c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
     plot!(plt_b, sol_i, vars=(1, 2), xlims=(-80.0, 0.0), c=:jet, line_z=test_rng[sol_idx], clims=(test_rng[1], test_rng[end]), legend=false)
end
plot(plt_a, plt_b, layout=(1, 2))

#%% Now lets run through several parameters using the wave model
nx = ny = 50
#b = Binomial(1, 0.75) #This means 75% of cells are sensitive to GABA
#null = Array{Float64}(rand(b, nx, ny))
param = :g_Ca
for val in LinRange(5.0, 10.0, 10) #This is the range of values
     print("[$(now())]: Setting up binomal nullification... ")
     println(" [$(now())]: Completed")
     net = GABA_PDE#(dU, U, p, t) -> GABA_PDE(dU, U, p, t, null)
     print("[$(now())]: Loading Parameters... ")
     conds_dict = read_JSON("params/GABA_conds.json")
     u0 = extract_dict(conds_dict, GABA_conds, dims=(nx, ny))
     pars_dict = read_JSON("params/GABA_params.json")
     #CHANGE ANY PARAMETER HERE
     pars_dict[param] = val
     p = extract_dict(pars_dict, GABA_pars)
     tspan = (0.0, 120e3)
     println(" [$(now())]: Completed")

     print("[$(now())]: Warming up model... ")
     prob = SDEProblem(net, noise, u0, tspan, p)
     @time warmup = solve(prob, SOSRI(), abstol=2e-2, reltol=2e-2, maxiters=1e7, save_everystep=false, progress=true, progress_steps=1)
     println(" [$(now())]: Completed")

     print("[$(now())]: Running model... ")
     prob = SDEProblem(net, noise, warmup[end], tspan, p)
     @time NetSol = solve(prob, SROCK1(), dt=1.0, abstol=2e-2, reltol=0.2, maxiters=1e7, progress=true, progress_steps=1, save_idxs=[1:(nx*ny)...])
     println(" [$(now())]: Completed")
     println(size(NetSol))
     # Step 7: Animate the solution
     print("[$(now())]: Animating solution... ")
     loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Modelling\GABA_experiments"
     animate_dt = 60.0
     anim = @animate for t = 1.0:animate_dt:NetSol.t[end]
          println("[$(now())]: Animating simulation...")
          frame_i = reshape(NetSol(t) |> Array, (nx, ny))
          heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
     end
     gif(anim, "$(loc)\\$(param)_$(round(val, digits = 1))animation.gif", fps=1000.0 / animate_dt)
     println(" [$(now())]: Completed")
end


#parameter spaces: 
#1) gGABA = [0.1 -> 1.0]
#2) gACh = [0.1 -> 0.215]
#2) gCa = [7.5 -> 8.0]
#3) gK = [1.0 -> 3.0]
#4) μ = [0.7-0.8]
