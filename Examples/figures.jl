using RetinalChaos
using Plots, Formatting, Colors
using LaTeXStrings
using Plots.Measures
import RetinalChaos: M_INF, N_INF, Φ, ħ
import RetinalChaos: diffuse
import RetinalChaos: calculate_threshold, max_interval_algorithim, timescale_analysis
#Define plotting attributes
font_title = Plots.font("Arial", 24)
font_axis = Plots.font("Arial", 12)
font_legend = Plots.font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
homedir = pwd()
pars_path = joinpath(homedir, "Settings")
figs_path = joinpath(homedir, "Examples")
#Set up default graph colors
v_color = :deepskyblue
n_color = :magenta
c_color = :green
a_color = :purple
b_color = :red
ach_color = :blue

#TODO: find out how to properly use xticks

#%% Setting up basic Settings
#Read JSON files for initial conditions and parameters
u0 = read_JSON(joinpath(pars_path,"conds.json")) |> extract_dict;
#Read JSON files for parameters and set the applied current to 10pA
p_dict = read_JSON(joinpath(pars_path, "params.json")) 
p_dict[:I_app] = 10.0
p = p_dict |> extract_dict
#Set the time span from 0s -> 60s and the sampling interval to 0.1
tspan = (0.0, 60e3); dt = 0.1 
#Create the problem
prob = ODEProblem(T_ode, u0, tspan, p)
@time sol = solve(prob); 

#%% Run all code here to generate Figure 1

###### Figure 1 A
xlims = (23.8e3, 24e3)
ticks = collect(xlims[1]:20:xlims[2])
labels = ["\$$(round(Int, x-xlims[1]))\$" for x in ticks]  
fig1_Aa = plot(sol, vars = [:v, :n,], plotdensity = 10000, tspan = xlims,
    c = [v_color n_color], lw = 2.0, 
    ylabel = ["\$v_t \$(mV)" "\$n_t\$"], xlabel = ["" "Time (ms)"], 
    xlims = xlims, xticks = (ticks, labels),
    layout = grid(2, 1)
)
    
vt_sol = map(t -> sol(t)[1], collect(xlims[1]:dt:xlims[2]))
nt_sol = map(t -> sol(t)[2], collect(xlims[1]:dt:xlims[2]))
fig1_Ab = plot(vt_sol, nt_sol, c = :black, 
    xlabel = "\$v_t\$ (mV)", ylabel = "\$n_t\$",
    title = "\$v_t\$ to \$n_t\$ Oscillation Diagram", titlefontsize = 12.0,
)

fig1_A = plot(fig1_Aa, fig1_Ab, 
    title = ["A" ""], title_location = :left, 
    layout = grid(1,2, widths = [0.7, 0.3])
)

##### Figure 1 B
v_rng = collect(-100:1.0:0)
f(v) = (p_dict[:C_0]+p_dict[:δ]*(-p_dict[:g_Ca] * M_INF(v, p_dict[:V1], p_dict[:V2]) * (v - p_dict[:E_Ca])) - (p_dict[:λ]*p_dict[:C_0]))/p_dict[:τc]
dc_rng = map(f, v_rng);


xlims = (22e3, 26e3) 
ticks = collect(xlims[1]:500:xlims[2])
labels = ["\$$(round(Int, x-xlims[1]))\$" for x in ticks]

fig1_Ba = plot(sol, vars = [:v, :c], tspan = xlims, plotdensity = 10000,
    c = [v_color c_color], lw = 2.0,
    ylabel = ["\$V_t\$ (mV)" "\$C_t\$ (mM)"], xlabel = ["" "Time (ms)"], xlims = xlims, 
    xticks = (ticks, labels),
    layout = grid(2, 1) 
)

fig1_Bb = plot(v_rng, dc_rng*1000, 
    c = c_color, linestyle = :dash, lw = 2.0, fill = (0, 0.5, c_color),
    xlabel = "\$V_t\$ (mV)",ylabel = "Inst. \$C_t\$ (mM/s)",
    title = "Instantaneous Calcium Change", titlefontsize = 12.0
)

fig1_B = plot(fig1_Ba, fig1_Bb, 
    title = ["B" "" "" ""], title_pos = :left,
    layout = grid(1, 2, widths=[0.7, 0.3])
)

###### Figure 1 C
c_rng = collect(0.0:0.01:1.5)
fa(c) = ((p_dict[:α]*c^4)/p_dict[:τa])
a_rng = map(fa, c_rng).*1000;

a_rng_2 = collect(0.0:0.01:1.5)
fb(a) = ((p_dict[:β]*a^4)/p_dict[:τb])
b_rng = map(fb, a_rng_2).*1000;

xlims = (0.0, 40e3)
ticks = collect(xlims[1]:5000:xlims[2])
labels = ["\$$(round(Int, (x-xlims[1])/1000))\$" for x in ticks]

fig1_Ca = plot(sol, vars = [:c, :a, :b, :v], 
    c = [c_color a_color b_color v_color], lw = 2.0, 
    ylabel = ["\$C_t\$ (mM)" "\$A_t\$ (mM)" "\$B_t\$ (mM)" "\$V_t\$ (mM)" ], xlabel = ["" "" "" "Time (s)"] , 
    xlims = xlims, xticks = (ticks, labels), ylims = [(0.0, 1.0) (0.0, 1.0) (0.0, 1.0) (-90.0, 0.0)],
    layout = grid(4,1)
)
fig1_Cb = plot(c_rng, a_rng, 
    c = a_color, lw = 2.0, linestyle = :dash, fill = (0, 0.2, a_color),
    xlabel = "\$C_t\$ (mM)", ylabel = "Inst. \$A_t\$", 
    xlims = (0.0, 0.6), ylims = (0.0, 1.0),
    title = "\$C_t\$ and \$A_t\$", titlefontsize = 10.0
)
fig1_Cc = plot(a_rng_2, b_rng,  
    c = b_color, lw = 2.0, linestyle = :dash, fill = (0, 0.2, b_color),
    xlabel = "\$A_t\$", ylabel = "Inst. \$B_t\$",
    xlims = (0.0, 0.6), ylims = (0.0, 1.0),
    title = "\$A_t\$ and \$B_t\$", titlefontsize = 10.0,    
)

fig1_C = plot(fig1_Ca, fig1_Cb, fig1_Cc, 
    title = ["C" "" "" "" "" ""], title_pos = :left,
    layout = grid(1,3, widths = [0.6, 0.2, 0.2])
)
###### Putting all parts together
fig1 = plot(fig1_A, fig1_B, fig1_C, 
    grid = false, legend = :none,
    layout = grid(3, 1, heights = [0.2, 0.3, 0.5]), 
    dpi = 300, size = (2250, 2625)
)

#Files need to be in .tiff or .eps format
savefig(fig1, joinpath(figs_path, "Figure1_ModelDynamics.svg"))
println("Figure 1 generated")

#%% Figure 2
#Extract the second burst timestamps for plotting
trace = map(t -> sol(t)[1], collect(tspan[1]:dt:tspan[end]))
sim_thresh = calculate_threshold(trace)
spike_array = (trace .> sim_thresh);
output, _, _, _ = max_interval_algorithim(spike_array; dt = dt);
burst2 = output[2]

v_rng = LinRange(-90.0, 10.0, 100)
fach(v) = Φ(v, p_dict[:k], p_dict[:V0])
ACh_rng = map(fach, v_rng);

xlims = (burst2[1] - 5e3, burst2[2]+5e3) 
xticks = (collect(xlims[1]:500:xlims[2]), collect(0:0.5:round(Int, tspan[end]/1000)))
frame_stops = LinRange(burst2[1], burst2[2], 4)
ach_stops = map(t->sol(t)[6], frame_stops)

fig2_Aa = plot(sol, vars = [:v, :e], c = [v_color ach_color], lw = 2.0,
    ylabel = ["\$V_t\$ (mV)" "\$E_t\$ (mM)"], xlabel = ["" "time (s)"], 
    xlims = xlims, 
    #xticks = xticks 
    legend = :none, layout = grid(2, 1), 
)

plot!(fig2_Aa[2], frame_stops, ach_stops, label = "Frame stops",
    seriestype = :scatter, marker = :star, markersize = 10.0)

fig2_Ab = plot(v_rng, ACh_rng, 
    c = ach_color, linestyle = :dot, lw = 4.0, 
    title = "Acetylcholine Release", titlefontsize = 12.0,
    xlabel = "\$V_t\$ (mV)",ylabel = "Inst. \$E_t\$ (mM/s)", label = "")

fig2_A = plot(fig2_Aa, fig2_Ab, 
    title = ["A" ""], titlepos = :left,
    layout = grid(1, 2, widths = [0.75, 0.25]), size = (700, 350))

#%%
#Create the PDE equations
nx, ny = (50, 50)
c1x, c1y = (round(Int, nx/2), round(Int, ny/4))
c2x, c2y = (round(Int, nx/2), round(Int, ny/4*3))
D = p_dict[:D]
lattice = zeros(nx, ny)
lattice[c1x, c1y] = 0.6
bp_model = Network(nx, ny)

time_range = collect(xlims[1]:xlims[2])
println(length(time_range))
lattice_c = zeros(nx, ny, length(time_range))
#Step through a simple diffusion
for (idx, t) in enumerate(time_range)
    #Each step, the cell releases 0.005 ACh
    if idx == 1
        lattice_c[c1x, c1y, idx] = sol(t)[6]
    else
        lattice_c[:, :, idx] = diffuse(lattice_c[:,:,idx-1], D, bp_model)
        lattice_c[c1x, c1y, idx] = sol(t)[6]
    end
end

fig2_B = plot(xlims = (0, nx), ylims = (0, ny), layout = grid(1, length(frame_stops)), grid = false, size = (1000, 250), 
    margin = 1mm
)

upper_lim = 0.08
for (idx, frame) in enumerate(frame_stops)
    frame_idx = round(Int, frame - xlims[1])
    #p = plot(xlims = (0.1, nx), ylims = (0.1, ny), 
    #    xaxis = false, yaxis = false)
    if idx == length(frame_stops)
        heatmap!(fig2_B[idx], lattice_c[:,:,frame_idx], c = :thermal, clims = (0.0, upper_lim), 
            aspect_ratio = :equal, grid = false,
            xaxis = false, yaxis = falseheat
        )
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :cyan)), label = "")
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :circle, c = :cyan, markersize = 10.0, label = "")
        scatter!(fig2_B[idx], [c2y], [c2x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :red)), label = "")
        scatter!(fig2_B[idx], [c2y], [c2x], marker = :circle, c = :red, markersize = 10.0, label = "")
        annotate!(fig2_B[idx], [40], [3], "t = $(time_range[frame_idx+1]/1000)s", :white)
    else
        heatmap!(fig2_B[idx], lattice_c[:,:,frame_idx], c = :thermal, clims = (0.0, upper_lim), 
            aspect_ratio = :equal, grid = false, colorbar = false, 
            xaxis = false, yaxis = false
        )
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :cyan)), label = "")
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :circle, c = :cyan, markersize = 10.0, label = "")
        scatter!(fig2_B[idx], [c2y], [c2x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :red)), label = "")
        scatter!(fig2_B[idx], [c2y], [c2x], marker = :circle, c = :red, markersize = 10.0, label = "")
        annotate!(fig2_B[idx], [40], [3], "t=$(time_range[frame_idx+1]/1000)s", :white)
    end
end
fig2_B
title!(fig2_B[1], "B", title_pos = :left)

fig2_Ca = plot(layout = (2,1))
plot!(fig2_Ca[1], sol, vars = [:e], c = :blue, label = "Cell 1 release",
    xlims = xlims, 
    #xticks = xticks, 
    ylabel = "\$E_t\$ (mM)",
    xlabel = ""
)

ach_rng = LinRange(0.001, 10e4, 100)
i_rng = map(a -> ħ(a, p_dict[:k_d]), ach_rng);
ach_external = lattice_c[c1x, c1y+1, :];
i_synaptic = map(a ->ħ(a, p_dict[:k_d]), ach_external);


plot!(fig2_Ca[1], time_range.+xlims[1], ach_external, label = "Cell 2 Extracellular", 
    c = :red, lw = 3.0, legend = :topleft)
ach_rel = map(t->sol(t)[6], frame_stops)
ach_ext = map(t->lattice_c[c1x, c1y+1, round(Int, t-xlims[1])], frame_stops)
plot!(fig2_Ca[1], frame_stops, ach_rel, label = "",
    seriestype = :scatter, marker = :star, markersize = 10.0)
plot!(fig2_Ca[1], frame_stops, ach_ext, label = "",
    seriestype = :scatter, marker = :star, markersize = 10.0)

plot!(fig2_Ca[2], time_range./1000.0, i_synaptic, label = "", 
    c = :purple, lw = 2.0,
    ylabel = "Norm. Current", xlabel = "time (s)", ylims = (0,1)
)


fig2_Cb = plot(ach_rng, i_rng, label = "",
    xlabel = "Extracellular [ACh] (mM)", ylabel = "Normalized Current", 
    title = "Induced \$I_{ACh}\$", titlefontsize = 12.0, 
    c = :green, lw = 3.0, 
)
fig2_C = plot(fig2_Ca, fig2_Cb, layout = grid(1,2, widths = [0.75, 0.25]))
title!(fig2_C[1], "C", titlepos = :left)

fig2 = plot(fig2_A, fig2_B, fig2_C, layout = grid(3, 1), size = (1000,1000))
#savefig(fig2, joinpath(figs_path, "Figure2_Acetylcholine_Rel_Diff_Activate.png"))