using Revise
using RetinalChaos #Load the package
using BSON, JLD2
#%%Lets figure out how to extract the waves
import RetinalChaos.param_path
import RetinalChaos.DiscreteCallback
#%% Load a model and make a Current injection step
using DiffEqCallbacks
tspan = (0.0, 10e3)
params_file = joinpath(param_path, "params.json")
conds_file = joinpath(param_path, "conds.json")
u0 = read_JSON(conds_file)|> extract_dict;
p = read_JSON(params_file) 
p[:g_ACh] = 0.0
p0 = p |> extract_dict

#Lets set up a callback for the current step
step_begin = 500.0
duration = 10.0
level = 1.0
cb = RetinalChaos.IC_callback(step_begin, duration, level)

#step_begin2 = step_begin + duration
#duration2 = 1000.0
#level2 = -100.0
#cb2 = RetinalChaos.IC_callback(step_begin2, duration2, level2; level_delta = 10.0)
#cb_set = RetinalChaos.CallbackSet(cb, cb2)

prob_func(i) = SDEProblem(T_sde, noise, u0 , tspan, p0, callback = cb);
prob = prob_func(10)
@time sol = solve(prob, SOSRI(), progress = true)
plot(sol, vars = 1)
#%% Set up a ensemble function
n_sims = 50
test_rng = repeat([p[:σ]], n_sims)
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, :σ, test_rng)
#make the ensemble problem and only save the voltage
simProb = EnsembleProblem(prob, prob_func = prob_func);
#%%
@time sim = solve(simProb, trajectories = n_sims, save_idxs = 1, EnsembleThreads(), callback = i -> cb_i(i));

#%%
plot(sim[1], vars = 1)
#%%
plot(sol, vars = 1)
plot!([step_begin], st = :vline)
plot!([step_begin+duration], st = :vline)
plot!([step_begin2], st = :vline)
plot!([step_begin2+duration2], st = :vline)


#%%
load_path = "F:\\Data\\Modelling\\mu_0\\"
p_dict = read_JSON("$(load_path)\\params.json", is_type = Dict{Symbol, Float32})
u_dict = read_JSON("$(load_path)\\iconds.json", is_type = Dict{Symbol, Float32})
sol = load_model(load_path, p_dict, u_dict, gpu = false)
#%%
RetinalChaos.save_solution(sol, load_path)
#%%
sol_loaded = RetinalChaos.load_solution(load_path)
