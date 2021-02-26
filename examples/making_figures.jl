#%% Running and analyzing the model using RetinalChaos.jl
using RetinalChaos
import RetinalChaos: Φ, diffuse
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

#%% Figure2 Acetylcholine Propagation dynamics
p = read_JSON("params.json") 
u0 = read_JSON("conds.json");
p[:I_app] = 2.0
p[:D] = 0.01
tspan = (0.0, 60e3);
prob = ODEProblem(T_ode, u0|> extract_dict, tspan, p|> extract_dict)
println("Time it took to simulate 60s:")
@time sol = solve(prob, abstol = 2e-2, reltol = 2e-2, maxiters = 1e7); 

plot(sol, vars = [:v, :e], layout = grid(2, 1))


v_rng = LinRange(-90.0, 10.0, 100);
AChi = 0.0
ρ = p[:ρ]
k = p[:k]
println(k)
V0 = p_dict[:V0]
τACh = p_dict[:τACh]
fach(v) = (ρ * Φ(v, k, V0)) - (AChi)/τACh
ACh_rng = map(fach, v_rng);

xlims = (38e3, 42e3) 
xticks = (collect(xlims[1]:500:xlims[2]), collect(0:0.5:round(Int, tspan[end]/1000)))

fig2_Aa = plot(sol, 
    vars = [:v, :e], 
    ylabel = ["\$V_t\$ (mV)" "\$E_t\$ (mM)"], 
    xlabel = ["" "time (s)"], 
    c = [v_color ach_color], lw = 2.0,
    legend = :none, 
    layout = grid(2, 1), 
    xlims = xlims, 
    xticks = xticks
)
frame_stops = [1.5e3, 2.0e3, 2.5e3, 3.0e3]
frame_stops .+= xlims[1]
ach_stops = map(t->sol(t)[6], frame_stops)
plot!(fig2_Aa[2], frame_stops, ach_stops, label = "Frame stops",
    seriestype = :scatter, marker = :star, markersize = 10.0)

fig2_Ab = plot(v_rng, ACh_rng, 
    c = ach_color, linestyle = :dot, lw = 4.0, 
    title = "Acetylcholine Release", titlefontsize = 12.0,
    xlabel = "\$V_t\$ (mV)",ylabel = "Inst. \$E_t\$ (mM/s)", label = "")

fig2_A = plot(fig2_Aa, fig2_Ab, 
    layout = grid(1, 2, widths = [0.75, 0.25]), size = (700, 350))
title!(fig2_A[1], "A", titlepos = :left)

nx, ny = (50, 50)
c1x, c1y = (round(Int, nx/2), round(Int, ny/4))
c2x, c2y = (round(Int, nx/2), round(Int, ny/4*3))
D = p_dict[:D]
lattice = zeros(nx, ny)
lattice[c1x, c1y] = 0.6
bp_model = Network(nx, ny)