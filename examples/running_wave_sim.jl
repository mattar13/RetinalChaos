using RetinalChaos

#Run a wave model 
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

println("File set up")
#%% Set up initial conditions
nx = 64; ny = 64
p = read_JSON(params_file);
p[:g_HCN] = 0.75
p[:g_ACh] = 0.50
p[:τa] = 15e3
p[:τb] = 15e3
p = extract_dict(p, tar_pars);
u0_network = extract_dict(read_JSON(conds_file), tar_conds, (nx, ny));
net = Network(nx, ny; μ = 0.40, version = :gHCN) #μ is the probability a cell is capable of being active
var = 1 #This is the variable we are interested in
save_idxs = [var*1:var*(nx*ny)...] #This is a list of indexes of all the variables we want to plot
#%% Lets warm up the solution first
PDEprob = SDEProblem(net, noise, u0_network, (0.0, 60e3), p)
print("Warming up solution:")
@time PDEsol = solve(PDEprob, SOSRA(), 
    abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
    save_everystep = false, progress = true, progress_steps = 1
    )
# extract the last solution and use it as the initial condition for future models 
# This may take a long time to run
print("Running model:")
ui_network = PDEsol[end]
PDEprob = SDEProblem(net, noise, ui_network, (0.0, 120e3), p)
@time PDEsol = solve(PDEprob, SOSRA(), 
    abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
    save_idxs = save_idxs, progress = true, progress_steps = 1
    )
#%%
plt = plot()
for i in 1:20
    println(i)
    plot!(plt, PDEsol.t, PDEsol[i,:])
end 
plt
#%% Plot as a .gif
t_rng = collect(1.0:50.0:PDEsol.t[end])
anim = @animate for t = t_rng
    println("Animating frame $t")
    frame_i = reshape(PDEsol(t), (nx, ny))
    heatmap(frame_i, ratio = :equal, grid = false,
            xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
            c = :curl, clims = (-70.0, 0.0),
    )
end
#%%
gif(anim, "examples\\gHCN_warmed_up2.gif", fps = 20)
#%%
contour(net.null)