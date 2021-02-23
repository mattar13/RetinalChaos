#%% Running and analyzing the model using RetinalChaos.jl
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

#%% Making Figure 1
#Set up parameters and other items
p_dict = read_JSON(params_file) 
p_dict[:I_app] = 10.0
p_dict[:g_ACh] = 0.0
p = p_dict |> extract_dict;
#Setup initial conditions
u0 = read_JSON(conds_file) |> extract_dict;
#Setup figure 1 problem
tspan = (0.0, 60e3)
prob = ODEProblem(T_ode, u0, tspan, p);
#Run the simulation
print("Time it took to simulate $(tspan[2]/1000)s:")
@time sol = solve(prob); 

#to do the analysis we should set a specific dt
dt = 0.1 #set the time differential
t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
sim_thresh = calculate_threshold(v_t)
spike_array = (v_t .> sim_thresh)
timestamps = get_timestamps(spike_array; dt = dt)
burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array, dt = dt);
ts_analysis = timescale_analysis(v_t, dt = dt)
spike_dur = sum(ts_analysis[1])/length(ts_analysis[1])

#%% Make figure 1
xlims = (burst_idxs[2][1], burst_idxs[2][1]+200)
elapsed_time = xlims[2]-xlims[1]
dx_lims = 20
xticks = (collect(xlims[1]:dx_lims:xlims[2]), Int64.(collect(0:dx_lims:elapsed_time)))
fig1_A = plot(sol, vars = [:v, :n,], legend = :none, plotdensity = Int(100e3),
    ylabel = ["\$V_t\$ (mV)" "\$N_t\$"], xlabel = ["" "time (ms)"], 
    lw = 2.0, c = [v_color n_color],
    layout = grid(2, 1),grid = false,
    xlims = xlims, xticks = xticks
)

xlims = (burst_idxs[2][1]-1500, burst_idxs[2][2]+1500)
elapsed_time = xlims[2] - xlims[1]
dx_lims = 500
xticks = (collect(xlims[1]:dx_lims:xlims[2]), (collect(0:dx_lims/1000:elapsed_time/1000)))
fig1_B = plot(sol, vars = [:v, :c], legend = :none, 
    ylabel = ["\$V_t\$ (mV)" "\$C_t\$"], xlabel = ["" "time (s)"], 
    lw = 2.0, c = [v_color c_color], 
    layout = grid(2, 1),grid = false,
    xlims = xlims, xticks = xticks
)

dx_lims = 5e3
xticks = (
    collect(sol.t[1]:dx_lims:sol.t[end]), 
    collect(sol.t[1]/1000:dx_lims/1000:sol.t[end]/1000)
    )
fig1_C = plot(sol, vars = [:c, :a, :b, :v], legend = :none, 
    ylabel = ["\$V_t\$ (mV)" "\$C_t\$"], xlabel = ["" "" "" "time (s)"], 
    lw = 2.0, c = [c_color a_color b_color v_color],
    layout = grid(4, 1),grid = false,
    xticks = xticks 
)

fig1 = plot(fig1_A, fig1_B, fig1_C, layout = grid(3,1, heights = [0.2, 0.3, 0.5]), size = (1000, 1000))

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
sfig1 = plot(dt_rng, spike_durs, label = "",
    xlabel = ["" "" "dt (ms)"], ylabel = ["Spike dur (ms)" "Burst dur (ms)" "IBI (ms)"], 
    xaxis = :log, 
    layout = (3,1))
plot!(sfig1[2], dt_rng, burst_durs, label = "")
plot!(sfig1[3], dt_rng, IBI_durs, label = "")
savefig(sfig1, joinpath(save_figs, "Supp_fig1.png"))