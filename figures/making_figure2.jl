#%% Running and analyzing the model using RetinalChaos.jl
using RetinalChaos
import RetinalChaos: Φ, diffuse, ħ
using Dates
using StatsBase, Statistics
using LaTeXStrings
#Setup the fonts and stuff

font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
#RetinalChaos.pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
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


#%% Figure2 Acetylcholine Propagation dynamics
p = read_JSON(params_file) 
u0 = read_JSON(conds_file);
p[:I_app] = 10.0
p[:D] = 0.01
tspan = (0.0, 60e3);
prob = ODEProblem(T_ode, u0|> extract_dict, tspan, p|> extract_dict)
println("Time it took to simulate 60s:")
@time sol = solve(prob, abstol = 2e-2, reltol = 2e-2, maxiters = 1e7); 
# Run an anlysis and extract the bursts
dt = 0.1 #set the time differential
t_rng = collect(tspan[1]:dt:tspan[2]) #set the time range
v_t = map(t -> sol(t)[1], t_rng); #extract according to the interval
sim_thresh = calculate_threshold(v_t)
spike_array = (v_t .> sim_thresh)
burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(spike_array, dt = dt);
burst_lims = burst_idxs[2]

dx_lims = 2000.0
xlims = (burst_lims[1]-2000, burst_lims[2]+2000)
elapsed_time = xlims[2]-xlims[1]
xticks = (collect(xlims[1]:dx_lims:xlims[2]), (collect(0:dx_lims/1000:elapsed_time/1000)))

#This is for drawing the diffusion frames
frame_stops = LinRange(burst_lims[1], burst_lims[2], 3)
ach_stops = map(t->sol(t)[6], frame_stops)
# Figure 2A
fig2_Aa = plot(sol, vars = [:v, :e], layout = grid(2, 1), plotdensity = Int64(200e3),
    c = [v_color e_color], legend = false, 
    xlabel = ["" "Time (s)"], ylabel = ["\$V_t\$ (mV)" "\$E_t\$ (nM)"], 
    xlims = xlims, xticks = xticks
)
#Plot the frame stops
plot!(fig2_Aa[2], frame_stops, ach_stops, label = "Frame stops",
    seriestype = :scatter, marker = :star, markersize = 10.0)
fig2_Ab = plot(v -> Φ(v, p[:k], p[:V0]), -90, 10, legend = false, xlabel = "\$V_t\$", ylabel = "\$\\Phi\$ (Normalized ACh release)")
fig2_A = plot(fig2_Aa, fig2_Ab, layout = grid(1,2))
title!(fig2_A[1], "A", titlepos = :left)

#Make the diffusion model
nx, ny = (50, 50)
c1x, c1y = (round(Int, nx/2), round(Int, ny/2))
c2x, c2y = (round(Int, nx/2), round(Int, ny/4*3))
D = p[:D] #this is the dX per dt
lattice = zeros(nx, ny)
lattice[c1x, c1y] = 0.6
#%%
bp_model = Network(nx, ny, μ = 0.15)
#%% Plot the binomial Nullification
binom = heatmap(bp_model.null, ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :grays, clims = (0.0, 1.0), 
)
title!(binom, "μ = 0.15")
savefig(binom, "figures\\binomial_15.png")
#%%
dt = 1.0
time_range = collect(xlims[1]:dt:xlims[2])
lattice_c = zeros(nx, ny, length(time_range))
for (idx, t) in enumerate(time_range)
    if idx == 1
        lattice_c[c1x, c1y, idx] = sol(t)[6]
    else
        lattice_c[:, :, idx] = diffuse(lattice_c[:,:,idx-1], D, bp_model)
        lattice_c[c1x, c1y, idx] = sol(t)[6]
    end
end

#%% This is a supplemental gif for what diffusion looks like
d_frame = 17
anim = @animate for frame = 1:d_frame:length(time_range)
    println("Animating frame $frame")
    frame_i = lattice_c[:,:, frame]
    heatmap(frame_i, ratio = :equal, grid = false,
            xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
            c = :thermal, clims = (0.0, 0.10), 
    )
    scatter!([c1y], [c1x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :cyan)), label = "")
    scatter!([c1y], [c1x], marker = :circle, c = :cyan, markersize = 10.0, label = "")
    #scatter!([c2y], [c2x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :red)), label = "")
    #scatter!([c2y], [c2x], marker = :circle, c = :red, markersize = 10.0, label = "")
    annotate!([40], [3], "t = $(frame/1000) ms", :white)
end
gif(anim, "figures\\ACh_Release.gif", fps = 1000/d_frame)
#%%
fig2_B = plot(xlims = (0, nx), ylims = (0, ny), layout = grid(1, length(frame_stops)), size = (1000, 250), 
    #margin = 1mm
)
upper_lim = 0.08
for (idx, frame) in enumerate(frame_stops)
    frame_idx = round(Int, frame-xlims[1])
    time = round((time_range[frame_idx+1] - frame_stops[1])/1000, digits = 1)
    println(time)
    if idx == length(frame_stops)
        heatmap!(fig2_B[idx], lattice_c[:,:,frame_idx], c = :thermal, clims = (0.0, upper_lim), 
            aspect_ratio = :equal, grid = false, colorbar_title = "\$E_t\$ (nM)",
            xaxis = false, yaxis = false
        )
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :cyan)), label = "")
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :circle, c = :cyan, markersize = 10.0, label = "")
        #scatter!(fig2_B[idx], [c2y], [c2x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :red)), label = "")
        #scatter!(fig2_B[idx], [c2y], [c2x], marker = :circle, c = :red, markersize = 10.0, label = "")
        annotate!(fig2_B[idx], [40], [3], "t = $(time) ms", :white)
    else
        heatmap!(fig2_B[idx], lattice_c[:,:,frame_idx], c = :thermal, clims = (0.0, upper_lim), 
            aspect_ratio = :equal, grid = false, colorbar = false, 
            xaxis = false, yaxis = false
        )
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :cyan)), label = "")
        scatter!(fig2_B[idx], [c1y], [c1x], marker = :circle, c = :cyan, markersize = 10.0, label = "")
        #scatter!(fig2_B[idx], [c2y], [c2x], marker = :hexagon, m = (125.0, :transparent, stroke(3.0, :red)), label = "")
        #scatter!(fig2_B[idx], [c2y], [c2x], marker = :circle, c = :red, markersize = 10.0, label = "")
        annotate!(fig2_B[idx], [40], [3], "t=$(time) ms", :white)
    end
end
title!(fig2_B[1], "B", title_pos = :left)

#%% Plotting Acetylcholine extracellular conc vs 

#%%
ach_external = lattice_c[c1x, c1y+1, :];
i_synaptic = map(a ->ħ(a, p[:k_d]), ach_external);

#%%
fig2_Ca = plot(layout = (2,1))
plot!(fig2_Ca[1], sol, vars = [:e], c = :blue, label = "Cell 1 release",
    xlims = xlims, 
    xticks = xticks, 
    ylabel = "\$E_t\$ (mM)",
    xlabel = ""
)
plot!(fig2_Ca[1], time_range.+xlims[1], ach_external, label = "Cell 2 Extracellular", 
    c = :red, lw = 3.0, legend = :topleft)
#%%
ach_rel = map(t->sol(t)[6], frame_stops)
ach_ext = map(t->lattice_c[c1x, c1y+1, Int(t-xlims[1])], frame_stops)
plot!(fig2_Ca[1], frame_stops, ach_rel, label = "",
    seriestype = :scatter, marker = :star, markersize = 10.0)
plot!(fig2_Ca[1], frame_stops, ach_ext, label = "",
    seriestype = :scatter, marker = :star, markersize = 10.0)

plot!(fig2_Ca[2], time_range./1000.0, i_synaptic, label = "", 
    c = :purple, lw = 2.0,
    ylabel = "Norm. Current", xlabel = "time (s)", ylims = (0,1)
)


fig2_Cb = plot(a -> ħ(a, p[:k_d]), 0.001, 6.0, xlabel = "\$E_t\$ (nM)", ylabel = "nAChR Activation (H)")

fig2_C = plot(fig2_Ca, fig2_Cb, layout = grid(1,2, widths = [0.75, 0.25]))
title!(fig2_C[1], "C", titlepos = :left)

fig2 = plot(fig2_A, fig2_B, fig2_C, layout = grid(3, 1), size = (1000,1000))