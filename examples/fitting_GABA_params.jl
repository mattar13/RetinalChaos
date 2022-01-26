using Revise
using RetinalChaos #Load the package
using BSON, JLD2
using LinearAlgebra
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

nx = 10 
ny = 100
DX = (0.2, 1.8) #Bias right and left
DY = (0.9, 1.1) #This breaks down to diffusion bias up and down

StencilE = RetinalChaos.diffusion_stencil(nx, ny; dX = DX, dY = DY)
a = zeros(nx, ny);
a[5, 50] = 100.0;
ai = RetinalChaos.diffuse(StencilE, a, 0.1)
aii = RetinalChaos.diffuse(StencilE, ai, 0.1)

p1 = Plots.heatmap(a, ratio = :equal);
p2 = Plots.heatmap(ai, ratio = :equal);
p3 = Plots.heatmap(aii, ratio = :equal);
plot(p1, p2, p3, layout = (1, 3))
