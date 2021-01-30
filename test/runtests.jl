#I will put all of the testing stuff here
using Revise #We put this here so even while we are editing we don't need to wait for the updates
using RetinalChaos

#%%
println("Testing opening parameters")
param_root = "params\\"
p = read_JSON(joinpath(param_root, "params.json")) |> extract_dict
u0 = read_JSON(joinpath(param_root, "conds.json")) |> extract_dict
tspan = (0.0, 1000.0)
println("Model components opened properly")

#%% Testing model components
println("Testing ODE")
ODEprob = ODEProblem(T_ode, u0, tspan, p);
@time ODEsol = solve(ODEprob, progress = true);
println("Ordinary Differential Equation working")

#%%
println("Testing SDE")
SDEprob = SDEProblem(T_sde, u0, tspan, p);
@time SDEsol = solve(SDEprob, SOSRI(), progress = true);
println("Stochastic Differential Equation working")