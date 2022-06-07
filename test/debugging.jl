using Revise
using RetinalChaos
using ABFReader

#We need to double check the properties of the physiological data




#%% lets try to make a save callback
results = RetinalChaos.SavedValues(Float64, Vector{Float64})
cb = RetinalChaos.SavingCallback(
    (u, t, integrator) -> reshape(u[:, :, 1], size(u, 1) * size(u, 2)),
    results
)

p_dict[:t_warm] = 300e3
p_dict[:t_run] = 1000.0
p0 = p_dict |> extract_dict
u0 = extract_dict(u_dict, p_dict[:nx], p_dict[:ny]) |> cu
gpu = true
version = :gACh
net = Network(p_dict[:nx], p_dict[:ny]; μ=p_dict[:μ], version=version, gpu=gpu)
NetProb = SDEProblem(net, noise, u0, (0.0f0, p_dict[:t_warm]), p0)
print("[$(now())]: Warming up the solution... ")

@time sol = solve(NetProb, SOSRI(),
    callback=cb,
    abstol=2e-2, reltol=2e-2, maxiters=1e7,
    progress=true, progress_steps=1,
    save_everystep=false
)
#We can the null out the solution because GPU is expensive
#%%
results.saveval
#first we want to save 

#%% Experiments with Current Clamp callbacks
step_begin = 500.0
duration = 10.0
level = 1.0
cb = RetinalChaos.IC_callback(step_begin, duration, level)

prob = SDEProblem(T_sde, noise, u0, tspan, p0, callback=cb);
@time sol = solve(prob, SOSRI(), progress=true)
plot(sol, vars=1)
#%% Set up a ensemble function
n_sims = 50
test_rng = repeat([p[:σ]], n_sims)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, :σ, test_rng)
#make the ensemble problem and only save the voltage
simProb = EnsembleProblem(prob, prob_func=prob_func);
#%%
@time sim = solve(simProb, trajectories=n_sims, save_idxs=1, EnsembleThreads(), callback=i -> cb_i(i));

#%%
plot(sim[1], vars=1)
#%%
plot(sol, vars=1)
plot!([step_begin], st=:vline)
plot!([step_begin + duration], st=:vline)
plot!([step_begin2], st=:vline)
plot!([step_begin2 + duration2], st=:vline)


#%%
load_path = "F:\\Data\\Modelling\\mu_0\\"
p_dict = read_JSON("$(load_path)\\params.json", is_type=Dict{Symbol,Float32})
u_dict = read_JSON("$(load_path)\\iconds.json", is_type=Dict{Symbol,Float32})
sol = load_model(load_path, p_dict, u_dict, gpu=false)
#%%
RetinalChaos.save_solution(sol, load_path)
#%%
sol_loaded = RetinalChaos.load_solution(load_path)
