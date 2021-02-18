#I will put all of the testing stuff here
using Revise #We put this here so even while we are editing we don't need to wait for the updates
using RetinalChaos

#%%
println("Testing opening parameters")
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")
p = read_JSON(params_file) |> extract_dict
u0_single = read_JSON(conds_file) |> extract_dict
tspan = (0.0, 1000.0)
println("Model components opened properly")

#%% Testing running Single cell ODE and SDEs
println("Testing ODE")
ODEprob = ODEProblem(T_ode, u0_single, tspan, p);
@time ODEsol = solve(ODEprob, progress = true);
println("Ordinary Differential Equation working")

#%% SDE
tspan = (0.0, 120e3)
println("Testing SDE")
SDEprob = SDEProblem(T_sde, u0_single, tspan, p);
@time SDEsol = solve(SDEprob, SOSRI(), progress = true, save_idxs = [1]);
println("Stochastic Differential Equation working")
#%%
using Plots
plot(SDEsol)
#%% Running a wave model
dx = 0.01; nx = 64; ny = 64;
#In a network all of the parameters and conditions need to be read differently
p_dict = read_JSON(params_file)
p = extract_dict(p_dict, tar_pars);

u_dict = read_JSON(conds_file);
u0_network = extract_dict(u_dict, tar_conds, (nx, ny));
#Establish the PDE network
net = Network(nx, ny; μ = 0.65)

var = 1 #This is the variable we are interested in
save_idxs = [var*1:var*(nx*ny)...] #This is a list of indexes of all the variables we want to plot

PDEprob = SDEProblem(net, noise, u0_network, (0.0, 6.0), p)
#%% This next step is quite lengthy, so before running this, make sure you have time
@time PDEsol = solve(PDEprob, SOSRI(), save_idxs = save_idxs, progress = true, progress_steps = 1)
#%% Animation
t_rng = collect(0.0:500.0:60e3)
anim = @animate for t = t_rng
    println("Animating frame $t")
    frame_i = reshape(PDEsol(t), (nx, ny))
    heatmap(frame_i, ratio = :equal, grid = false,
            xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
            c = :curl, clims = (-70.0, 0.0),
    )
end
#%%
gif(anim, "plot.gif", fps = 20)