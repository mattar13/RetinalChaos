#Using this will set the default plotting font
#font_title = Plots.font("Arial", 24)
#font_axis = Plots.font("Arial", 12)
#font_legend = Plots.font("Arial", 8)
#plotlyjs(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
import Plots.@animate

export equilibria_object
const v_color = :deepskyblue
const n_color = :magenta
const c_color = :green
const a_color = :purple
const b_color = :red
const e_color = :blue
const w_color = :gray
export v_color, n_color, c_color, a_color, b_color, e_color, w_color


#=function animate_solution(sol, save_path::String;
    animate_dt=60.0, verbose=true
)
    nx = ny = Int64(sqrt(size(sol, 1)))
    #we can use this to build a solution without GPU
    println("[$(now())]: Animating simulation...")
    anim = @animate for t = 1.0:animate_dt:sol.t[end]

        frame_i = reshape(sol(t) |> Array, (nx, ny))
        heatmap(frame_i, ratio=:equal, grid=false,
            xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny),
            c=:curl, clims=(-70.0, 0.0),
        )
    end
    gif(anim, "$(save_path)\\animation.gif", fps=1000.0 / animate_dt)
end
=#

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

@recipe function f(c1::codim_object{1,T}; vars=:v, scatter=false) where {T<:Real}
    #var_idx = findall(x -> x==vars, tar_conds)[1]
    var_idx = vars |> u_find
    points = map(x -> x[1], c1.points)
    sorted_dims = sortperm(points)
    points = points[sorted_dims]
    saddle_p = map(eq -> length(eq.saddle) > 0 ? eq.saddle[1][var_idx] : NaN, c1.equilibria)[sorted_dims]
    stable_p = map(eq -> length(eq.stable) > 0 ? eq.stable[1][var_idx] : NaN, c1.equilibria)[sorted_dims]
    unstable_p = map(eq -> length(eq.unstable) > 0 ? eq.unstable[1][var_idx] : NaN, c1.equilibria)[sorted_dims]
    unstable_focus_p = map(eq -> length(eq.unstable_focus) > 0 ? eq.unstable_focus[1][var_idx] : NaN, c1.equilibria)[sorted_dims]
    stable_focus_p = map(eq -> length(eq.stable_focus) > 0 ? eq.stable_focus[1][var_idx] : NaN, c1.equilibria)[sorted_dims]

    plotted_stable = false
    plotted_unstable = false
    plotted_saddle = false
    plotted_unstable_focus = false
    plotted_stable_focus = false

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
