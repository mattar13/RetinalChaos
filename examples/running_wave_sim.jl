using RetinalChaos

#Run a wave model 
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

#%% Set up initial conditions
nx = 96; ny = 96
p = read_JSON(params_file)
p[:g_Ca] = 10.0
p[:g_K] = 5.0
p = extract_dict(p, tar_pars);
u0_network = extract_dict(read_JSON(conds_file), tar_conds, (nx, ny));
nx = 96; ny = 96
net = Network(nx, ny; μ = 0.15, version = :Sparse_Activity) #μ is the probability a cell is capable of being active
var = 1 #This is the variable we are interested in
save_idxs = [var*1:var*(nx*ny)...] #This is a list of indexes of all the variables we want to plot
PDEprob = SDEProblem(net, noise, u0_network, (0.0, 60e3), p)
#%% This may take a long time to run
@time PDEsol = solve(PDEprob, SOSRI(), save_idxs = save_idxs, progress = true, progress_steps = 1)

#%% Plot as a .gif
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
gif(anim, "examples\\plot.gif", fps = 20)
