using Revise
using RetinalChaos
using ABFReader

#Lets try to fix the PDE stuff so that we can adjust certain variables

#%% Take info from Zhous 2004 paper
#extract the function for ACh activation
import RetinalChaos.ħ
p = read_JSON(params_file)
#x -> voltage
#p -> [gACh, k_d]
model1(x, par) = map(v -> par[1] * ħ(1.0, par[2]) * (v - p[:E_ACh]), x)
v_known = [35.0, 15.0, -5.0, -25.0, -45.0, -75.0] #Zhou known v's
i_known = [0.0, 0.0, -25.0, -50.0, -100.0, -200.0] #Zhou known I's
p0 = [0.75, 0.1]
fit = curve_fit(model, v_known, i_known, p0, lower=[0.0, 0.0], upper=[Inf, Inf])
p1 = plot(v -> model1(v, fit.param), -75, 35)
p2 = plot(x -> ħ(x, fit.param[2]), 0.001, 0.1)
plot(p1, p2, layout=grid(2, 1))
#%% Lets measure the difference between the adaptive properties of the model vs the timestepping 
p = read_JSON(params_file)
#%%
p[:I_app] = 10.0
p[:g_ACh] = 0.0
u0 = read_JSON(conds_file);
#Setup figure problem
tspan = (0.0, 60e3)
prob_sde = SDEProblem(T_sde, u0 |> extract_dict, tspan, p |> extract_dict);
print("Time it took to simulate $(tspan[2]/1000)s:")
@time sol = solve(prob_sde);

dt = 0.1
timesteps = tspan[1]:dt:tspan[2]
timestep_average = sum(sol(timesteps, idxs=1)) / length(timesteps)
adaptive_average = sum(sol(sol.t, idxs=1)) / length(sol.t)
plot(sol, vars=[:v])
hline!([timestep_average], c=:black, label="timestep average")
hline!([adaptive_average], c=:green, label="adaptive average")

#%% Lets run a loop to test what the best and fastest threshold is
import RetinalChaos.calculate_threshold
times = Float64[]
results = Float64[]
dts = LinRange(100.0, 500.0, 200)
for dt in dts
    println("Running test for threshold dt: $dt")
    thresh = calculate_threshold(model, dt=dt)
    time = @elapsed calculate_threshold(model, dt=dt)
    push!(times, time)
    push!(results, thresh[1])
end
#%%
p1 = plot(dts, results, st=:scatter)
p2 = plot(dts, times, st=:scatter)
plot(p1, p2, layout=grid(2, 1))


#%% This supplemental figure compares the dt to the analysis accuracy
dt_rng = range(0.005, 0.10, length=50)
spike_durs = Float64[];
spike_dur_stds = Float64[];
burst_durs = Float64[];
burst_dur_stds = Float64[];
IBI_durs = Float64[];
IBI_dur_stds = Float64[];

for dt in dt_rng
    println("testing dt= $dt")
    t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
    v_t = map(t -> sol(t)[1], t_rng) #extract according to the interval
    print("Analysis took:")
    @time ts_analysis = timescale_analysis(v_t, dt=dt)

    spike_dur = sum(ts_analysis[1]) / length(ts_analysis[1])
    spike_dur_std = std(ts_analysis[1])
    push!(spike_durs, spike_dur)
    push!(spike_dur_stds, spike_dur_std)

    burst_dur = sum(ts_analysis[2]) / length(ts_analysis[2])
    burst_dur_std = std(ts_analysis[2])
    push!(burst_durs, burst_dur)
    push!(burst_dur_stds, burst_dur_std)

    IBI_dur = sum(ts_analysis[3]) / length(ts_analysis[3])
    IBI_dur_std = std(ts_analysis[3])
    push!(IBI_durs, IBI_dur)
    push!(IBI_dur_stds, IBI_dur_std)
end
sfig2 = plot(dt_rng, spike_durs, label="",
    xlabel=["" "" "dt (ms)"], ylabel=["Spike dur (ms)" "Burst dur (ms)" "IBI (ms)"],
    xaxis=:log,
    layout=(3, 1))
plot!(sfig1[2], dt_rng, burst_durs, label="")
plot!(sfig1[3], dt_rng, IBI_durs, label="")
savefig(sfig1, joinpath(save_figs, "Supp_fig1.png"))

#%% Supplemental figure. Analysis of parameter change
test_rng = range(0.5, 2.0, length=100) #this ranges from halving the parameter to doubling it
p = read_JSON(params_file); #Because of my catch, we can keep these as dictionaries 
u0 = read_JSON(conds_file);
#We can set the initial voltage sensitivty here 
p[:σ] = 10.0
tspan = (0.0, 60e3)
prob = ODEProblem(T_ode, u0, tspan, p);
results = zeros(length(T_ode.ps), length(test_rng), 3)
#Walk through each parameter
for (par_idx, par) in enumerate(T_ode.ps)
    println("Simulating test range for parameters: $par")
    #Setup an ensemble function to test each parameter
    par_rng = test_rng * p_dict[Symbol(par)]
    prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, par_rng)
    ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

    print("Time it took to simulate $(tspan[2]/1000)s:")
    @time sim = solve(ensemble_prob, EnsembleThreads(), trajectories=length(test_rng))
    println("Simulation for parameter $par completed")
    for (sol_idx, sol) in enumerate(sim)
        #println(sol |> length)
        dt = 0.1 #set the time differential according to supp figure 1
        t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
        v_t = map(t -> sol(t)[1], t_rng) #extract according to the interval
        ts_analysis = timescale_analysis(v_t, dt=dt)
        for i = 1:length(ts_analysis)
            results[par_idx, sol_idx, i] = sum(ts_analysis[i]) / length(ts_analysis[i])
        end
    end
    println("Statistics for $par calculated")
end
#%%
for (idx, p) in enumerate(T_ode.ps)
    println("Plotting results for $p")
    sfig3i = plot(layout=grid(3, 1), size=(1000, 800))
    plot!(sfig3i[1], test_rng, results[idx, :, 1], label="$p", xlabel="Norm Par", ylabel="Spike Duration")
    plot!(sfig3i[2], test_rng, results[idx, :, 2] ./ 1000, label="$p", xlabel="Norm Par", ylabel="Burst Duration")
    plot!(sfig3i[3], test_rng, results[idx, :, 3] ./ 1000, label="$p", xlabel="Norm Par", ylabel="IBI")
    savefig(sfig3i, joinpath(save_figs, "gradient_for_param_$p.png"))
end
#sfig3i



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
