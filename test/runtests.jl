#I will put all of the testing stuff here
using Revise #We put this here so even while we are editing we don't need to wait for the updates
using RetinalChaos

#%%
println("Testing opening parameters")
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")
p = read_JSON(params_file)
u0_single = read_JSON(conds_file)
tspan = (0.0, 60e3)
println("Model components opened properly")
#%% Testing running Single cell ODE and SDEs
println("Testing ODE")
ODEprob = ODEProblem(T_ode, u0_single|>extract_dict, tspan, p|>extract_dict);
print("ODE run takes:")
@time ODEsol = solve(ODEprob, progress = true);
println("Ordinary Differential Equation working")
#%%
plot!(ODEsol, vars = [:v])
#%% SDE
println("Testing SDE")
SDEprob = SDEProblem(T_sde, u0_single|>extract_dict, tspan, p|>extract_dict);
print("SDE run takes:")
@time SDEsol = solve(SDEprob, SOSRI(), progress = true, save_idxs = [1]);
println("Stochastic Differential Equation working")
plot(SDEsol);
println("Plotting Solution successful")
#%% Running a wave model
dx = 0.01; nx = 64; ny = 64;
#In a network all of the parameters and conditions need to be read differently

p_dict = read_JSON(params_file); #Open the parameters as a dictionary
p = extract_dict(p_dict, tar_pars); #Extract the parameters into a usuable form

u_dict = read_JSON(conds_file); #Open the conditions as a dictionary
u0_network = extract_dict(u_dict, tar_conds, (nx, ny)); #Extract the conditions into a usuable form

#Establish the PDE network
net = Network(nx, ny; μ = 0.65)
var = 1 #This is the variable we are interested in
save_idxs = [var*1:var*(nx*ny)...] #This is a list of indexes of all the variables we want to plot
#Make the problem
PDEprob = SDEProblem(net, noise, u0_network, tspan, p)

#%% This next step is quite lengthy, so before running this, make sure you have time
print("SDE/PDE run takes:")
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
println("Animation successful")

#%% Testing Area
#%% Running a comparison between HCN+ cell and HCN- cell
p = read_JSON(params_file)
p[:g_ACh] = 0.0
u0 = read_JSON(conds_file); tspan = (0.0, 120e3)
HCNp_prob = SDEProblem(T_sde, u0|>extract_dict, tspan, p_HCNp|>extract_dict);
@time HCNp_sol = solve(HCNp_prob, SOSRI(), abstol = 2e-2, reltol = 2e-2, maxiters = 1e7, progress = true);

#%%
mA =
bA 
f(A) = 