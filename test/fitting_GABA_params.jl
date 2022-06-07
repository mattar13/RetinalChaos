using Revise
using RetinalChaos #Load the package
using LinearAlgebra
using StaticArrays
#%%Lets figure out how to extract the waves
import RetinalChaos: param_path, GABA_conds, GABA_pars, GABA_ODE, GABA_SDE

params_file = joinpath(param_path, "GABA_params.json")
conds_file = joinpath(param_path, "GABA_conds.json")
#%% Running a GABA 1D model
p_dict = read_JSON(params_file)
u_dict = read_JSON(conds_file)
#p_dict[:g_GABA] = 2.0
p = extract_dict(p_dict, GABA_pars) 
u0 = extract_dict(u_dict, GABA_conds)
tspan = (0.0, 300e3)
probODE = ODEProblem(GABA_ODE, u0, tspan, p);
probSDE = SDEProblem(GABA_SDE, noise, u0, tspan, p);


@time solODE = solve(probODE, progress = true);
@time solSDE = solve(probSDE, progress = true);
#%% Plotting results
plt_a = plot(solSDE, vars = [1, 6, 7], layout = (3, 1))

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
     abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
     save_idxs = [1:(nx*ny)...],
     progress = true, progress_steps = 1
)

#%%
dFrame = 100.0
anim = @animate for t in NetSol.t[1]:dFrame:NetSol.t[end]
     println(t)
     frame = reshape(NetSol(t), (nx, ny))
     plot(frame, st = :heatmap, ratio = :equal, c = :jet, clims = (-70.0, 10.0))
end

gif(anim, fps = 10.0)

plot(Array(NetSol)[2,:])