#Because the wave interface can get kinda slow in Julia, we have to run it in atom
using DifferentialEquations, Sundials
using Plots, Colors, LaTeXStrings
using ProgressLogging
include("RetinalChaos.jl");
cd("C:\\Users\\mtarc\\JuliaScripts\\RetinalChaos\\Notebooks")

#Defining Initial conditions and parameters
p_dict = read_JSON("params.json", is_type = Dict{Symbol,Float64});
u_dict = read_JSON("conds.json", is_type = Dict{Symbol,Float64});
p_dict[:g_ACh] = 0.0
tspan = (0.0, 60e3)
#Run 3 seperate trials
#Trial 1
sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
append_modeldata("DataSheet.xlsx", df_stats, df_params)
#Trial 2
sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
append_modeldata("DataSheet.xlsx", df_stats, df_params)
#Trial 3
sol_array, df_params, df_stats = run_model(p_dict, u_dict, tspan)
append_modeldata("DataSheet.xlsx", df_stats, df_params)

wave_arr = SDE_sol_arr[:, :, 1, :]
threshold = calculate_threshold(wave_arr)
spike_arr = wave_arr .>= threshold
θr = find_maxISI(spike_arr; dt = 10.0)
burst_arr = convolve_bursts(spike_arr, θr; dt = 10.0)
nx, ny, t = size(wave_arr)
spike_raster = reshape(spike_arr, (nx * ny, t))
burst_raster = convolve_bursts(spike_raster, θr; dt = 10.0)

#run the model
raster = reshape(SDE_sol_arr, (nx * ny, 7, size(SDE_sol_arr, 4)))
pick_ten = rand((1:size(raster, 1)), 10)
p = plot(layout = grid(6, 1), cbar = :none, legend = false);
plot!(p[1], raster[pick_ten, 1, :]', lw = 4.0, c = :blues, line_z = pick_ten');
plot!(p[2], raster[pick_ten, 2, :]', lw = 4.0, c = :PuRd, line_z = pick_ten');
plot!(p[3], raster[pick_ten, 3, :]', lw = 4.0, c = :greens, line_z = pick_ten', clims = pick_ten');
plot!(p[4], 1 .- raster[pick_ten, 4, :]', lw = 4.0, c = :BuPu, line_z = pick_ten', clims = pick_ten');
plot!(p[5], raster[pick_ten, 5, :]', lw = 4.0, c = :reds, line_z = pick_ten', clims = pick_ten');
plot!(p[6], raster[pick_ten, 6, :]', lw = 4.0, c = :bgy, line_z = pick_ten', clims = pick_ten')

p = plot(layout = grid(7, 1), size = (2000, 2000));
heatmap!(p[1], raster[:, 1, :], c = :curl);
heatmap!(p[2], burst_raster, c = :grayscale);
heatmap!(p[3], raster[:, 2, :], c = :RdPu);
heatmap!(p[4], raster[:, 3, :], c = :kgy);
heatmap!(p[5], 1 .- raster[:, 4, :], c = :BuPu);
heatmap!(p[6], raster[:, 5, :], c = :reds);
heatmap!(p[7], raster[:, 6, :], c = :bgy)

#Graphing trials
anim = @animate for i = 1:2:size(SDE_sol_arr, 4)
    println("Animating frame $i")
    heatmap(SDE_sol_arr[:, :, 1, i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :curl, clims = (-70.0, 0.0),
    )
end
gif(anim, "vt_map.gif", fps = 50)
mp4(anim, "vt_map.mp4", fps = 50)

anim = @animate for i = 1:2:size(SDE_sol_arr, 4)
    println("Animating frame $i")
    heatmap(burst_arr[:,:,i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :grays, clims = (0.0, 1.0),
    )
    #contourf!(SDE_sol_arr[:,:,1,i])
end
gif(anim, "spike_map.gif", fps = 50.0)

anim = @animate for i = 1:size(SDE_sol_arr, 4)
    println("Animating frame $i")
    heatmap(
        SDE_sol_arr[:, :, 3, i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :kgy, clims = (0.0, 1.0),
    )
end
gif(anim, "Ca_map.gif", fps = 10.0)

m_ach = maximum(SDE_sol_arr[:, :, 6, :])
anim = @animate for i = 1:size(SDE_sol_arr, 4)
    println("Animating frame $i")
    heatmap( SDE_sol_arr[:, :, 6, i],ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :bgy, clims = (0.0,m_ach)
    )
end
gif(anim, "ACh_map.gif", fps = 60)
mp4(anim, "ACh_map.mp4", fps = 100)

m_ach = maximum(SDE_sol_arr[:, :, 6, :])
anim = @animate for i = 1:2:size(SDE_sol_arr, 4)
    println("Animating frame $i")
    p = plot(layout = grid(2, 2), size = (800, 800))
    heatmap!(p[1], SDE_sol_arr[:, :, 1, i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :curl, clims = (-70.0, 0.0),
    )
    heatmap!(p[2], burst_arr[:,:,i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :grays, clims = (0.0, 1.0),
    )
    heatmap!(p[3], SDE_sol_arr[:, :, 3, i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :kgy, clims = (0.0, 1.0),
    )
    heatmap!(p[4], SDE_sol_arr[:, :, 6, i], ratio = :equal, grid = false,
        xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
        c = :bgy, clims = (0.0, m_ach),
    )
end
gif(anim, "FullMap.gif", fps = 50)
mp4(anim, "map.mp4", fps = 50)

keys = Tuple(BurstModel.params)
p_nt = NamedTuple{keys}(p0);
#write_JSON(p_nt, "params.json")
