#I will put all of the testing stuff here
using Revise
using RetinalChaos

#%%
println("Testing opening parameters")
param_root = "params\\"
p = read_JSON(joinpath(param_root, "params.json")) |> extract_dict
u0 = read_JSON(joinpath(param_root, "conds.json")) |> extract_dict
tspan = (0.0, 1000.0)
println("Model components opened properly")
#%% Testing model components
ODEprob = ODEProblem(T_ode, u0, tspan, p);
@time ODEsol = solve(ODEprob, progress = true);
println("Ordinary Differential Equation working")
#%% For some reason I can't understand, they swicth from ODE system symbols and variables to this
u0map = map(x -> u0[x], T_sde.states)
pmap = map(x -> x[2], p)
#%%
SDEprob = SDEProblem(RetinalChaos.T_sde, u0, tspan, p);
@time SDEsol = solve(SDEprob, SOSRI(), progress = true);
println("Stochastic Differential Equation working")