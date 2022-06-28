#%% Lets measure the difference between the adaptive properties of the model vs the timestepping 
pars_dict = read_JSON(params_file)
pars_dict[:I_app] = 10.0
pars_dict[:g_ACh] = 0.0
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
