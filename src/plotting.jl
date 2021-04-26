#Using this will set the default plotting font
#font_title = Plots.font("Arial", 24)
#font_axis = Plots.font("Arial", 12)
#font_legend = Plots.font("Arial", 8)
#plotlyjs(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)


const v_color = :deepskyblue
const n_color = :magenta
const c_color = :green
const a_color = :purple
const b_color = :red
const e_color = :blue
const w_color = :gray
export v_color, n_color, c_color, a_color, b_color, e_color, w_color

# ALL OF THESE SHOULD BE RECIPES SO WE DON'T HAVE TO IMPORT PLOTS WITH THE FXN
#function frame_draw(sol_array; idx = :all, saveas = :gif)
#    threshold = calculate_threshold(sol_array[:,:,1,:])
#    spike_arr = sol_array[:,:,1,:] .>= threshold
#    burst_arr = extract_burstmap(spike_array);
#    
#    nx, ny, var, t = size(sol_array)
#    if idx == :all
#        save_name = :Full
#        anim = @animate for i = 1:2:size(sol_array, 4)
#            println("Animating frame $i")
#            p = plot(layout = grid(2, 2), size = (800, 800))
#            heatmap!(p[1], sol_array[:, :, 1, i], ratio = :equal, grid = false,
#                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
#                c = :curl, clims = (-70.0, 0.0),
#            )
#            heatmap!(p[2], burst_arr[:,:,i], ratio = :equal, grid = false,
#                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
#                c = :grays, clims = (0.0, 1.0),
#            )
#            heatmap!(p[3], sol_array[:, :, 3, i], ratio = :equal, grid = false,
#                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
#                c = :kgy, clims = (0.0, 1.0),
#            )
#            heatmap!(p[4], sol_array[:, :, 6, i], ratio = :equal, grid = false,
#                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
#                c = :bgy, clims = (0.0, m_ach),
#            )
#        end
#    elseif idx == :spike
#        save_name = :spike
#        anim = @animate for i = 1:2:size(sol_array, 4)
#            println("Animating frame $i")
#            heatmap(burst_arr[:,:,i], ratio = :equal, grid = false,
#                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
#                c = :grays, clims = (0.0, 1.0),
#            )
#            #contourf!(SDE_sol_arr[:,:,1,i])
#        end
#    else
#        save_name = model_conds[idx]
#        anim = @animate for i = 1:2:size(sol_array, 4)
#            println("Animating frame $i")
#            heatmap(sol_array[:, :, idx, i], ratio = :equal, grid = false,
#                xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
#                c = :curl, clims = (minimum(sol_arr[:, :, idx, i]), maximum(sol_arr[:, :, idx, i]))
#            )
#        end
#    end
#
#    if saveas == :mp4
#        mp4(anim, "$(save_name)t_map.mp4", fps = 50)
#    elseif saveas == :gif
#        gif(anim, "$(save_name)t_map.gif", fps = 50)
#    end
#end

#function raster_plot(sol_array)
#    threshold = calculate_threshold(sol_array[:,:,1,:])
#    spike_arr = sol_array[:,:,1,:] .>= threshold
#    nx, ny, t = size(spike_arr)
#    spike_raster = reshape(spike_arr, (nx * ny, t))
#    burst_raster = convolve_bursts(spike_raster; dt = 10.0)
#    raster = reshape(sol_array, (nx * ny, 7, t))
#    p = plot(layout = grid(7, 1), size = (2000, 2000));
#    heatmap!(p[1], raster[:, 1, :], c = :curl);
#    heatmap!(p[2], burst_raster, c = :grayscale);
#    heatmap!(p[3], raster[:, 2, :], c = :RdPu);
#    heatmap!(p[4], raster[:, 3, :], c = :kgy);
#    heatmap!(p[5], 1 .- raster[:, 4, :], c = :BuPu);
#    heatmap!(p[6], raster[:, 5, :], c = :reds);
#    heatmap!(p[7], raster[:, 6, :], c = :bgy)
#    p
#end


#function trace_plot(sol_array)
#    nx, ny, var, t = size(sol_array)
#    raster = reshape(sol_array, (nx * ny, 7, size(sol_array, 4)))
#    pick_ten = rand((1:size(raster, 1)), 10)
#    p = plot(layout = grid(6, 1), cbar = :none, legend = false);
#    plot!(p[1], raster[pick_ten, 1, :]', lw = 4.0, c = :blues, line_z = pick_ten');
#    plot!(p[2], raster[pick_ten, 2, :]', lw = 4.0, c = :PuRd, line_z = pick_ten');
#    plot!(p[3], raster[pick_ten, 3, :]', lw = 4.0, c = :greens, line_z = pick_ten', clims = pick_ten');
#    plot!(p[4], 1 .- raster[pick_ten, 4, :]', lw = 4.0, c = :BuPu, line_z = pick_ten', clims = pick_ten');
#    plot!(p[5], raster[pick_ten, 5, :]', lw = 4.0, c = :reds, line_z = pick_ten', clims = pick_ten');
#    plot!(p[6], raster[pick_ten, 6, :]', lw = 4.0, c = :bgy, line_z = pick_ten', clims = pick_ten')
#    p
#end


@recipe function f(eq::equilibria_object; vars = [:v, :n])
    seriestype := :scatter
    markersize := 8.0
    if eq.stable != []
        @series begin
            label := "Stable"
            seriescolor := :green
            marker := :star
            x = []
            y = []
            for pt in eq.stable
                push!(x, pt[vars[1]|>u_find])
                push!(y, pt[vars[2]|>u_find])
            end
            x, y
        end
    end
    if eq.unstable != []
        @series begin
            label := "Unstable"
            seriescolor := :red
            marker := :star
            x = []
            y = []
            for pt in eq.unstable
                push!(x, pt[vars[1]|>u_find])
                push!(y, pt[vars[2]|>u_find])
            end
            x, y
        end
    end
    if eq.saddle != []
        @series begin
            label := "Saddle"
            seriescolor := :blue
            markershape = :star
            x = []
            y = []
            for pt in eq.saddle
                push!(x, pt[vars[1]|>u_find])
                push!(y, pt[vars[2]|>u_find])
            end
            x, y
        end
    end
    if eq.unstable_focus != []
        @series begin
            label := "Unstable Focus"
            seriescolor := :red
            x = []
            y = []
            for pt in eq.unstable_focus
                push!(x, pt[vars[1]|>u_find])
                push!(y, pt[vars[2]|>u_find])
            end
            x, y
        end
    end
    if eq.stable_focus != []
        @series begin
            label := "Stable Focus"
            seriescolor := :green
            x = []
            y = []
            for pt in eq.stable_focus
                push!(x, pt[vars[1]|>u_find])
                push!(y, pt[vars[2]|>u_find])
            end
            x, y
        end
    end
end

@recipe function f(z::Float64, eq::equilibria_object; 
        view = :xyz, legend = :none)
    legend := legend
    seriestype := :scatter
    markersize := 8.0
    if eq.stable != []
        @series begin
            label := ""
            seriescolor := :green
            marker := :star

            xs = []
            ys = []
            zs = []
            for pt in eq.stable
                push!(zs, z)
                push!(xs, pt[1])
                push!(ys, pt[2])
            end
           if view == :xyz
                xs, ys, zs
            elseif view == :xzy
                xs, zs, ys
            elseif view == :zxy
                zs, xs, ys
            elseif view == :zyx
                zs, ys, xs
            elseif view == :yzx
                ys, zs, xs
            elseif view == :yxz
                ys, xs, zs
            elseif view == :xy
                xs, ys
            elseif view == :xz
                xs, zs
            elseif view == :yz
                ys, zs
            elseif view == :yx
                ys, xs
            elseif view == :zx
                zs, xs
            else
                zs, ys
            end
        end
    end
    if eq.unstable != []
        @series begin
            label := ""
            seriescolor := :red
            marker := :star
            xs = []
            ys = []
            zs = []
            for pt in eq.unstable
                push!(zs, z)
                push!(xs, pt[1])
                push!(ys, pt[2])
            end
           if view == :xyz
                xs, ys, zs
            elseif view == :xzy
                xs, zs, ys
            elseif view == :zxy
                zs, xs, ys
            elseif view == :zyx
                zs, ys, xs
            elseif view == :yzx
                ys, zs, xs
            elseif view == :yxz
                ys, xs, zs
            elseif view == :xy
                xs, ys
            elseif view == :xz
                xs, zs
            elseif view == :yz
                ys, zs
            elseif view == :yx
                ys, xs
            elseif view == :zx
                zs, xs
            else
                zs, ys
            end
        end
    end
    if eq.saddle != []
        @series begin
            label := ""
            seriescolor := :blue
            markershape := :star
            xs = []
            ys = []
            zs = []
            for pt in eq.saddle
                push!(zs, z)
                push!(xs, pt[1])
                push!(ys, pt[2])
            end
            if view == :xyz
                xs, ys, zs
            elseif view == :xzy
                xs, zs, ys
            elseif view == :zxy
                zs, xs, ys
            elseif view == :zyx
                zs, ys, xs
            elseif view == :yzx
                ys, zs, xs
            elseif view == :yxz
                ys, xs, zs
            elseif view == :xy
                xs, ys
            elseif view == :xz
                xs, zs
            elseif view == :yz
                ys, zs
            elseif view == :yx
                ys, xs
            elseif view == :zx
                zs, xs
            else
                zs, ys
            end
        end
    end
    if eq.unstable_focus != []
        @series begin
            label := ""
            seriescolor := :red
            xs = []
            ys = []
            zs = []
            for pt in eq.unstable_focus
                push!(zs, z)
                push!(xs, pt[1])
                push!(ys, pt[2])
            end
            if view == :xyz
                xs, ys, zs
            elseif view == :xzy
                xs, zs, ys
            elseif view == :zxy
                zs, xs, ys
            elseif view == :zyx
                zs, ys, xs
            elseif view == :yzx
                ys, zs, xs
            elseif view == :yxz
                ys, xs, zs
            elseif view == :xy
                xs, ys
            elseif view == :xz
                xs, zs
            elseif view == :yz
                ys, zs
            elseif view == :yx
                ys, xs
            elseif view == :zx
                zs, xs
            else
                zs, ys
            end
        end
    end
    if eq.stable_focus != []
        @series begin
            label := "Stable Focus"
            seriescolor := :green
            xs = []
            ys = []
            zs = []
            for pt in eq.stable_focus
                push!(zs, z)
                push!(xs, pt[1])
                push!(ys, pt[2])
            end
            if view == :xyz
                xs, ys, zs
            elseif view == :xzy
                xs, zs, ys
            elseif view == :zxy
                zs, xs, ys
            elseif view == :zyx
                zs, ys, xs
            elseif view == :yzx
                ys, zs, xs
            elseif view == :yxz
                ys, xs, zs
            elseif view == :xy
                xs, ys
            elseif view == :xz
                xs, zs
            elseif view == :yz
                ys, zs
            elseif view == :yx
                ys, xs
            elseif view == :zx
                zs, xs
            else
                zs, ys
            end
        end
    end
end

function extract_equilibria(c2::codim_object{2, T}, eq_type::Symbol; eq_var::Int64 = 1, view = :xyz) where T <: Real
    xs = []; ys = []; zs = [];
    for idx = 1:length(c2.points)
        eq = c2.equilibria[idx]
        if isa(c2.vars, Symbol)
            x = c2.points[idx]
            y = 0.0
        else
            x,y = c2.points[idx]
        end
        if eq_type == :stable && eq.stable != []
            for eq_points in eq.stable
                push!(xs, x)
                push!(ys, y)
                push!(zs, eq_points[eq_var])
            end
        end
        if eq_type == :unstable && eq.unstable != []
            for eq_points in eq.unstable
                push!(xs, x)
                push!(ys, y)
                push!(zs, eq_points[eq_var])
            end
        end
        if eq_type == :saddle && eq.saddle != []
            for eq_points in eq.saddle
                push!(xs, x)
                push!(ys, y)
                push!(zs, eq_points[eq_var])
            end
        end
        if eq_type == :stable_focus && eq.stable_focus != []
            for eq_points in eq.stable_focus
                push!(xs, x)
                push!(ys, y)
                push!(zs, eq_points[eq_var])
            end
        end
        if eq_type == :unstable_focus && eq.unstable_focus != []
            for eq_points in eq.unstable_focus
                push!(xs, x)
                push!(ys, y)
                push!(zs, eq_points[eq_var])
            end
        end
    end 
    if view == :xyz
        xs, ys, zs
    elseif view == :xzy
        xs, zs, ys
    elseif view == :zxy
        zs, xs, ys
    elseif view == :zyx
        zs, ys, xs
    elseif view == :yzx
        ys, zs, xs
    elseif view == :yxz
        ys, xs, zs
    elseif view == :xy
        xs, ys
    elseif view == :xz
        xs, zs
    elseif view == :yz
        ys, zs
    elseif view == :yx
        ys, xs
    elseif view == :zx
        zs, xs
    else
        zs, ys
    end
end

@recipe function f(c1::codim_object{1, Float64}; vars = :v, scatter = false)
    #var_idx = findall(x -> x==vars, tar_conds)[1]
    var_idx = vars |> u_find
    points = map(x -> x[1], c1.points);
    
    saddle_p = map(eq -> length(eq.saddle) > 0 ? eq.saddle[1][var_idx] : NaN, c1.equilibria);
    stable_p = map(eq -> length(eq.stable) > 0 ? eq.stable[1][var_idx] : NaN, c1.equilibria);
    unstable_p = map(eq -> length(eq.unstable) > 0 ? eq.unstable[1][var_idx] : NaN, c1.equilibria);
    unstable_focus_p = map(eq -> length(eq.unstable_focus) > 0 ? eq.unstable_focus[1][var_idx] : NaN, c1.equilibria);
    stable_focus_p = map(eq -> length(eq.stable_focus) > 0 ? eq.stable_focus[1][var_idx] : NaN, c1.equilibria);
    
    plotted_stable = false; plotted_unstable = false; plotted_saddle = false
    plotted_unstable_focus = false; plotted_stable_focus = false

    if !all(isnan.(stable_p))
        @series begin
            if plotted_stable == false
                label := "Stable"
                plotted_stable = true
            else
                label := ""
            end     
            seriescolor := :green
            if scatter
                marker := :star
            end
            points, stable_p
        end
    end
    
    if !all(isnan.(unstable_p))
        @series begin
            if plotted_unstable == false
                label := "Unstable"
                plotted_unstable = true
            else
                label := ""
            end     
            seriescolor := :red
            if scatter
                marker := :star
            end
            points, unstable_p
        end
    end
    
    if !all(isnan.(saddle_p))
        @series begin
            if plotted_saddle == false
                label := "Saddle"
                plotted_saddle = true
            else
                label := ""
            end   
            seriescolor := :blue
            if scatter
                marker := :star
            end
            points, saddle_p
        end
    end
    
    if !all(isnan.(stable_focus_p))
        @series begin
            if plotted_stable_focus == false 
                label := "Stable Focus"
                plotted_stable_focus = true
            else
                label := ""
            end      
            linestyle := :dash
            seriescolor := :green
            if scatter
                marker := :circle
            end
            points, stable_focus_p
        end
    end
    
    if !all(isnan.(unstable_focus_p))
        @series begin
            if plotted_unstable_focus == false 
                label := "Unstable Focus"
                plotted_unstable_focus = true
            else
                label := ""
            end     
            linestyle := :dash
            seriescolor := :red
            if scatter
                marker := :circle
            end
            points, unstable_focus_p
        end
    end
end

@recipe function f(c2::codim_object{2,T}; eq_var = 1, view = :xyz) where T <: Real
    println("This one is called for some reason")
    #legend := :none
    markersize := 8.0
    #Plotting begins here
    if isa(c2.vars, Symbol)
        if view == :thresholds
            @series begin
                legend := true
                label := "Threshold"
                seriescolor := :blue
                lw := 4.0
                extract_thresholds(c1_map)
            end
        elseif view == :baselines
            @series begin
                legend := true
                label := "Baselines"
                seriescolor := :green
                lw := 4.0
                extract_baselines(c1_map)
            end
        else
            @series begin
                legend := true
                seriestype := :scatter
                label := "Stable"
                seriescolor := :green
                marker := :star
                extract_equilibria(c2, :stable; eq_var = eq_var, view = :xz)
            end

            @series begin
                seriestype := :scatter
                label := "Unstable"
                seriescolor := :red
                marker := :star
                extract_equilibria(c2, :unstable; eq_var = eq_var, view = :xz)
            end


            @series begin
                seriestype := :scatter
                label := "Saddle"
                seriescolor := :blue
                marker := :star
                extract_equilibria(c2, :saddle; eq_var = eq_var, view = :xz)
            end


            @series begin
                seriestype := :scatter
                label := "Stable_focus"
                seriescolor := :green
                marker := :circle
                extract_equilibria(c2, :stable_focus; eq_var = eq_var, view = :xz)
            end


            @series begin
                
                seriestype := :scatter
                label := "Unstable_focus"
                seriescolor := :red
                marker := :circle
                extract_equilibria(c2, :unstable_focus; eq_var = eq_var, view = :xz)
            end  
        end
    else
        if view == :thresholds

            xspan, yspan, threshmap = extract_thresholds(c2)
            null_map = (threshmap .!= Inf) .* Inf
            legend := true
            @series begin
                seriestype := :heatmap
                seriescolor := :thermal
                xspan, yspan, threshmap
            end
            @series begin
                seriestype := :heatmap
                seriescolor := :black
                xspan, yspan, null_map
            end
        elseif view == :baselines
            xspan, yspan, basemap = extract_baselines(c2)
            null_map = (basemap .!= Inf) .* Inf
            legend := true
            @series begin
                seriestype := :heatmap
                seriescolor := :thermal
                xspan, yspan, basemap
            end
            @series begin
                seriestype := :heatmap
                seriescolor := :black
                xspan, yspan, null_map
            end
        else    
            @series begin
                seriestype := :scatter
                label := "Stable"
                seriescolor := :green
                marker := :star
                extract_equilibria(c2, :stable;eq_var = eq_var, view = view)
            end

            @series begin
                seriestype := :scatter
                label := "Unstable"
                seriescolor := :red
                marker := :star
                extract_equilibria(c2, :unstable; eq_var = eq_var, view = view)
            end

            @series begin
                seriestype := :scatter
                label := "Saddle"
                seriescolor := :blue
                marker := :star
                extract_equilibria(c2, :saddle; eq_var = eq_var, view = view)
            end

            @series begin
                seriestype := :scatter
                label := "Stable_focus"
                seriescolor := :green
                marker := :circle
                extract_equilibria(c2, :stable_focus; eq_var = eq_var, view = view)
            end

            @series begin
                seriestype := :scatter
                label := "Unstable_focus"
                seriescolor := :red
                marker := :circle
                extract_equilibria(c2, :unstable_focus; eq_var = eq_var, view = view)
            end  
        end
    end
end
