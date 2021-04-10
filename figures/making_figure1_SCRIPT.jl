#%% Making the figure for the overlay of the model
using RetinalChaos #This might be all that we need for figure 1
#Setup the fonts and stuff
import RetinalChaos.M_INF

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
u0 = read_JSON(conds_file);
#Setup figure problem
tspan = (0.0, 120e3)
prob = ODEProblem(T_ode, u0 |> extract_dict, tspan, p |> extract_dict);

print("Time it took to simulate $(tspan[2]/1000)s:")
@time sol = solve(prob); 

#%% to do the analysis we should set a specific dt
dt = 0.1 #set the time differential
v_thresh = calculate_threshold(sol, dt = dt)
timestamps = get_timestamps(sol, dt = dt)
burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(sol, dt = dt);
ts_analysis = timescale_analysis(sol, dt = dt)
#%%
spike_dur = sum(ts_analysis[1])/length(ts_analysis[1])

#%% Figure part A
burst_lims = burst_idxs[2]
xlims = (burst_lims[1], burst_lims[1]+200)
elapsed_time = xlims[2]-xlims[1]
dx_lims = 20
xticks = (collect(xlims[1]:dx_lims:xlims[2]), Int64.(collect(0:dx_lims:elapsed_time)))


fig1_Aa = plot(sol, vars = [:v, :n,], label = "", plotdensity = Int(1e5),
    ylabel = ["Vₜ (mV)" "Nₜ"], xlabel = ["" "Time (ms)"], 
    lw = 2.0, c = [v_color n_color],
    layout = grid(2, 1),grid = false,
    xlims = xlims, xticks = xticks, 
    link = :x, margin = 0.0mm, widen = false
)
plot!(fig1_Aa[1], xticks = false)
hline!(fig1_Aa, [v_thresh], label = "Spike threshold", c = :red, linewidth = 2.0, legend = :bottomright)

t_rng = burst_idxs[2][1]:dt:burst_idxs[2][2]
fig1_Ab = plot(sol(t_rng, idxs = 2), sol(t_rng, idxs = 1), xlabel = "Nₜ", ylabel = "Vₜ (mV)", label = "", c = :black)
hline!(fig1_Ab, [v_thresh], label = "Spike threshold", seriestype = :scatter, c = :red, legend = :bottomright)

fig1_A = plot(fig1_Aa, fig1_Ab, layout = grid(1,2, widths = (0.75, 0.25)), margin = 1mm)
title!(fig1_A[1], "A", titlepos = :left)

# Figure 1 B
xlims = (burst_lims[1]-1500, burst_lims[2]+1500)
elapsed_time = xlims[2] - xlims[1]
dx_lims = 500
xticks = (collect(xlims[1]:dx_lims:xlims[2]), (collect(0:dx_lims/1000:elapsed_time/1000)))
fig1_Ba = plot(sol, vars = [:v, :c], legend = :none, plotdensity = Int64(1e5),
    ylabel = ["Vₜ (mV)" "[Cₜ] (μM)"], xlabel = ["" "Time (s)"], 
    lw = 2.0, c = [v_color c_color], 
    layout = grid(2, 1),grid = false,
    xlims = xlims, xticks = xticks, 
    margin = 0.0mm# added margins call here
)
plot!(fig1_Ba[1], xticks = false)
f(v) = p[:δ]*(-p[:g_Ca] * M_INF(v, p[:V1], p[:V2]) * (v - p[:E_Ca]))
fig1_Bb = plot(f, -50, 10.0, c = :green, label = "", linewidth = 3.0, linestyle = :dash, xlabel = "Vₜ (mV)", ylabel = "[δC] (mM)")
fig1_B = plot(fig1_Ba, fig1_Bb, layout = grid(1,2, widths = (0.75, 0.25)))

title!(fig1_B[1], "B", titlepos = :left)

# Figure 1 C
xlims = (burst_idxs[2][1], burst_idxs[4][1])
dx_lims = 10e3
xticks = (
    collect(xlims[1]:dx_lims:xlims[2]), 
    round.(Int64, collect((xlims[1]-xlims[1])/1000:dx_lims/1000:(xlims[2]-xlims[1])/1000))
    )
fig1_Ca = plot(sol, vars = [:c, :a, :b, :v], legend = :none, 
    ylabel = ["[Cₜ](μM)" "Aₜ" "Bₜ" "Vₜ (mV)"], xlabel = ["" "" "" "Time (s)"], 
    lw = 2.0, c = [c_color a_color b_color v_color],
    layout = grid(4, 1), grid = false,
    xlims = xlims,
    xticks = xticks, margin = 0.0mm 
)

plot!(fig1_Ca[1], xticks = false)
plot!(fig1_Ca[2], xticks = false)
plot!(fig1_Ca[3], xticks = false)

t_rng = xlims[1]:dt:xlims[2]
fig1_Cb = plot(sol(t_rng, idxs = 4), sol(t_rng, idxs = 3), label = "", ylabel = "[Cₜ] (μM)", xlabel = "Aₜ", c = a_color, linewidth = 3.0)
fig1_Cc = plot(sol(t_rng, idxs = 5), sol(t_rng, idxs = 4), label = "", ylabel = "Aₜ", xlabel = "Bₜ", c = b_color, linewidth = 3.0)
fig1_Cd = plot(sol(t_rng, idxs = 1), sol(t_rng, idxs = 5), label = "", ylabel = "Bₜ", xlabel = "Vₜ (mV)", c = v_color, linewidth = 3.0)
fig1_C = plot(fig1_Ca, plot(fig1_Cb, fig1_Cc, fig1_Cd, layout = grid(3,1)), layout = grid(1,2, widths = (0.75, 0.25)))

title!(fig1_C[1], "C", titlepos = :left)

fig1 = plot(fig1_A, fig1_B, fig1_C, 
    layout = grid(3,1, heights = [0.2, 0.3, 0.5]), 
    size = (1000, 1000), grid = false
    )
#%% Save the figure if you are satusfyed with it
savefig(fig1, joinpath(save_figs, "Fig1_Model_Dynamics.png"))