#%% Running and analyzing the model using RetinalChaos.jl
using RetinalChaos #This might be all that we need for figure 1
#using StatsBase, Statistics
#Setup the fonts and stuff

font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
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

#%% Adjust plotting settings
i_app = 10.0


#%% Making Figure 1
#Set up parameters and other items
p = read_JSON(params_file) 
p[:I_app] = 10.0
p[:g_ACh] = 0.0
#Setup initial conditions
u0 = read_JSON(conds_file);
#Setup figure 1 problem
tspan = (0.0, 60e3)
prob = ODEProblem(T_ode, u0 |> extract_dict, tspan, p |> extract_dict);
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

#%% Figure part A
xlims = (burst_idxs[2][1], burst_idxs[2][1]+200)
elapsed_time = xlims[2]-xlims[1]
dx_lims = 20
xticks = (collect(xlims[1]:dx_lims:xlims[2]), Int64.(collect(0:dx_lims:elapsed_time)))
fig1_Aa = plot(sol, vars = [:v, :n,], legend = :none, plotdensity = Int(100e3),
    ylabel = ["Vₜ (mV)" "Nₜ"], xlabel = ["" "Time (ms)"], 
    lw = 2.0, c = [v_color n_color],
    layout = grid(2, 1),grid = false,
    xlims = xlims, xticks = xticks
)

dV = 0.1
vrng = sol(burst_idxs[2][1]:dV:burst_idxs[2][2], idxs = 1)
dN = 0.1
nrng = sol(burst_idxs[2][1]:dN:burst_idxs[2][2], idxs = 2)
fig1_Ab = plot(vrng, nrng, ylabel = "Nₜ", xlabel = "Vₜ (mV)")
fig1_A = plot(fig1_Aa, fig1_Ab, layout = grid(1,2))
title!(fig1_A[1], "A", titlepos = :left)


#%% Figure part B
xlims = (burst_idxs[2][1]-1500, burst_idxs[2][2]+1500)
elapsed_time = xlims[2] - xlims[1]
dx_lims = 500
xticks = (collect(xlims[1]:dx_lims:xlims[2]), (collect(0:dx_lims/1000:elapsed_time/1000)))
fig1_B = plot(sol, vars = [:v, :c], legend = :none, plotdensity = Int64(100e3),
    ylabel = ["Vₜ (mV)" "[Cₜ] mM"], xlabel = ["" "Time (s)"], 
    lw = 2.0, c = [v_color c_color], 
    layout = grid(2, 1),grid = false,
    xlims = xlims, xticks = xticks
)
title!(fig1_B[1], "B", titlepos = :left)

dx_lims = 5e3
xticks = (
    collect(sol.t[1]:dx_lims:sol.t[end]), 
    collect(sol.t[1]/1000:dx_lims/1000:sol.t[end]/1000)
    )
fig1_C = plot(sol, vars = [:c, :a, :b, :v], legend = :none, 
    ylabel = ["Cₜ" "Aₜ" "Bₜ" "Vₜ"], xlabel = ["" "" "" "time (s)"], 
    lw = 2.0, c = [c_color a_color b_color v_color],
    layout = grid(4, 1),grid = false,
    xticks = xticks 
)
title!(fig1_C[1], "C", titlepos = :left)

fig1 = plot(fig1_A, fig1_B, fig1_C, layout = grid(3,1, heights = [0.2, 0.3, 0.5]), size = (1000, 1000), grid = false)
#%% Save the figure if you are satusfyed with it
savefig(fig1, joinpath(save_figs, "Fig1_Model_Dynamics.png"))