# This file contains some seperate data analysis and experiments I am running.

using RetinalChaos
using Dates
using StatsBase, Statistics
#Setup the fonts and stuff

font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)

#Set up the file root and default parameters
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

#save everything in the figures folder
save_figs = "figures\\"
if isdir(save_figs) == false
    #The directory does not exist, we have to make it 
    mkdir(save_figs)
end

#%% Running a noise single trace stepping through parameters
test_rng = range(0.0, 0.5, length = 100) #this ranges from halving the parameter to doubling it
p = read_JSON(params_file); #Because of my catch, we can keep these as dictionaries 
u0 = read_JSON(conds_file);

tspan = (0.0, 60e3)
prob = SDEProblem(T_sde, u0|>extract_dict, tspan, p|>extract_dict);
#print("Time it took to simulate $(tspan[2]/1000)s:")
#@time sol = solve(prob, SOSRI(), abstol = 2e-2, reltol = 2e-2, maxiters = 1e7); 
#plot(sol, vars = [:v, :W], layout = grid(2,1))
#%%
par_idx = findall(x -> x==:σ, Symbol.(T_sde.ps))
prob_func(prob, i, repeat) = ensemble_func(prob, i, repeat, par_idx, test_rng)
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func);
print("Running a ensemble simulation for :")
@time sim = solve(ensemble_prob, SOSRI(), abstol = 0.2, reltol = 2e-2, maxiters = 1e7, trajectories = length(test_rng), EnsembleThreads());
#%%
results = zeros(3, length(test_rng))
trace_plot = plot(legend = false)
for (sol_idx, sol) in enumerate(sim)
    plot!(trace_plot, sol, vars = [:v, :W], line_z = sol_idx, c = :jet, zlim = (1, length(sim)), layout = grid(2,1))
    #println(sol |> length)
    dt = 0.1 #set the time differential according to supp figure 1
    t_rng = collect(sol.t[1]:dt:sol.t[end]) #set the time range
    v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
    ts_analysis = timescale_analysis(v_t, dt = dt)
    for i = 1:length(ts_analysis)
        results[i, sol_idx] = sum(ts_analysis[i])/length(ts_analysis[i])
    end
end
trace_plot
#%%
p1 = plot(test_rng, results[1, :])
p2 = plot(test_rng, results[2, :])
p3 = plot(test_rng, results[3, :])
sfig1 = plot(p1, p2, p3, layout = grid(3,1))
sfig1
#%% This supplemental figure compares the dt to the analysis accuracy
dt_rng = range(0.005, 0.10, length = 50)
spike_durs = Float64[]; spike_dur_stds = Float64[]
burst_durs = Float64[]; burst_dur_stds = Float64[]
IBI_durs = Float64[]; IBI_dur_stds = Float64[]

for dt in dt_rng
    println("testing dt= $dt")
    t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
    v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
    print("Analysis took:")
    @time ts_analysis = timescale_analysis(v_t, dt = dt)

    spike_dur = sum(ts_analysis[1])/length(ts_analysis[1])
    spike_dur_std = std(ts_analysis[1])
    push!(spike_durs, spike_dur)
    push!(spike_dur_stds, spike_dur_std)
    
    burst_dur = sum(ts_analysis[2])/length(ts_analysis[2])
    burst_dur_std = std(ts_analysis[2])
    push!(burst_durs, burst_dur)
    push!(burst_dur_stds, burst_dur_std)

    IBI_dur = sum(ts_analysis[3])/length(ts_analysis[3])
    IBI_dur_std = std(ts_analysis[3])  
    push!(IBI_durs, IBI_dur)
    push!(IBI_dur_stds, IBI_dur_std)
end
sfig2 = plot(dt_rng, spike_durs, label = "",
    xlabel = ["" "" "dt (ms)"], ylabel = ["Spike dur (ms)" "Burst dur (ms)" "IBI (ms)"], 
    xaxis = :log, 
    layout = (3,1))
plot!(sfig1[2], dt_rng, burst_durs, label = "")
plot!(sfig1[3], dt_rng, IBI_durs, label = "")
savefig(sfig1, joinpath(save_figs, "Supp_fig1.png"))

#%% 