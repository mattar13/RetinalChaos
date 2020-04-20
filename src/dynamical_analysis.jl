#This function iteratively assigns variables to a parameter or conditions within a range
function ensemble_func(prob, i, repeat; pars = [:I_app], par_range = LinRange(0.5, 15.0, 100))
    new_pars = prob.p
    new_conds = prob.u0
    stochastic = false
    for var in pars
        conds = var |> u_find
        pars = var |> p_find

        if prob.p[(:σ |> p_find)[1]] > 0.0
            stochastic = true
        end
        if length(conds) > 0
            new_conds[conds[1]] = par_range[i]
        elseif length(pars) > 0
            new_pars[pars[1]] = par_range[i]
        end
    end
    if stochastic
        return SDEProblem(prob.f, noise, new_conds, prob.tspan, new_pars)
    else
        return ODEProblem(prob.f, new_conds, prob.tspan, new_pars)
    end
end


function phase_plane(ui, pi; vars = [:v, :n], xlims = (-90.0, 10.0), ylims = (-0.10, 5.0), resolution = 100)
    du = similar(ui)
    n_vars = length(ui)
    var_idx = [(vars[1]|>u_find)[1], (vars[2]|>u_find)[1]]
    phase_plane = zeros(resolution, resolution, n_vars)
    for (idx_x, x) in enumerate(LinRange(xlims[1], xlims[2], resolution))
        for (idx_y, y) in enumerate(LinRange(ylims[1], ylims[2], resolution))
            ui[var_idx].= (x, y)
            BurstModel(du, ui, pi, 0.0)
            for idx_var = 1:n_vars
                phase_plane[idx_x, idx_y, idx_var] = du[idx_var]
            end
        end
    end
    phase_plane
end

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

#First we make a clojure function of the model
function all_in_one(ui, pi)
    du = similar(ui)
    BurstModel(du, ui, pi, 0.0)
    du
end

struct equilibria_object{T}
    stable::Array{Array{T}}
    unstable::Array{Array{T}}
    saddle::Array{Array{T}}
    unstable_focus::Array{Array{T}}
    stable_focus::Array{Array{T}}
end

#Conduct a stability analysis of the current
function find_equilibria(ui, pi;
        vars = [:v, :n], xlims = (-90.0, 10.0), ylims = (-1.0, 5.0), resolution = 10,
        precision = 2, check_min = 1e-8
    )
    stable = Array{Array{Float64}}([])
    unstable = Array{Array{Float64}}([])
    saddle = Array{Array{Float64}}([])
    unstable_focus = Array{Array{Float64}}([])
    stable_focus = Array{Array{Float64}}([])
    storage = Array{Array{Float64}}([])

    var_idx = [(vars[1]|>u_find)[1], (vars[2]|>u_find)[1]]
    copy_u = copy(ui)
    for (idx_x, x) in enumerate(LinRange(xlims[1], xlims[2], resolution))
        #Iterate through the x range
        for (idx_y, y) in enumerate(LinRange(ylims[1], ylims[2], resolution))
            #Iterate through the y range looking for stable points
            #We really may only need to check for stable points if the dU is low
            copy_u[var_idx] .= (x, y)
            df(x) = all_in_one(x, pi)
            res = nlsolve(df, copy_u)
            #println(res)
            equilibria = map(x -> round(x, digits = precision), res.zero)
            check = df(res.zero)
            if equilibria in storage
                nothing
            elseif any(isnan, res.zero)
                nothing
            elseif check[1] > check_min || check[2] > check_min
                nothing
            else
                push!(storage, equilibria)
                #push!(eq_dict[:all], res.zero)
                #If the equilibria is not in the stable or unstable points
                #Test it's stability
                J_mat = ForwardDiff.jacobian(x->all_in_one(x, pi), res.zero)[[var_idx...], [var_idx...]]
                ev = eigvals(J_mat)

                if sign(real(ev[1])) != sign(real(ev[2]))
                    #If the sign of the eigenvalues are opposite, then it is a unstable saddle
                    push!(saddle, res.zero)
                else
                    if imag(ev[1]) != 0 #This detects whether or not the equilibria is a focus
                        if real(ev[1]) > 0.0
                            push!(unstable_focus, res.zero)
                        else
                            push!(stable_focus, res.zero)
                        end
                    else
                        if ev[1] > 0.0
                            push!(unstable, res.zero)
                        else
                            push!(stable, res.zero)
                        end
                    end
                end
            end
        end
    end
    return equilibria_object{Float64}(stable, unstable, saddle, unstable_focus, stable_focus)
end

struct codim_object{N, T}
    vars::NTuple{N, Symbol}
    points::Array{NTuple{N, T}}
    equilibria::Array{equilibria_object{T}}
end


"""
Make a codim object
"""
function codim_map(ui, pi, codim::Symbol;
        c1_lims = (-10.0, 10.0), resolution = 50, eq_res = 3,
    )
    uv = copy(ui);
    pv = copy(pi)
    points_list = Array{Tuple{Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], resolution)
    for (idx1, c1) in enumerate(c1_range)
        pv = copy(p0)
        uv = copy(u0)
        pv[(codim |> p_find)...] = c1
        equilibria = find_equilibria(uv, pv; resolution = eq_res)
        points = c1
        #println(points |> typeof)
        push!(points_list, (points,))
        push!(equilibria_list, equilibria)
    end
    codim_object((codim,), points_list, equilibria_list)
end

function codim_map(ui, pi, codim::Tuple{Symbol, Symbol};
        c1_lims = (-10.0, 10.0), c2_lims = (0.0, 10.0), resolution = 50, eq_res = 3,
    )
    uv = copy(ui);
    pv = copy(pi)
    points_list = Array{Tuple{Float64, Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], resolution)
    c2_range = LinRange(c2_lims[1], c2_lims[2], resolution)
    for (idx1, c1) in enumerate(c1_range)
        #println(idx1)
        for (idx2, c2) in enumerate(c2_range)
            pv = copy(p0)
            uv = copy(u0)
            pv[(codim[1] |> p_find)...] = c1
            pv[(codim[2] |> p_find)...] = c2
            equilibria = find_equilibria(uv, pv; resolution = eq_res)
            points = (c1, c2)
            #println(points |> typeof)
            push!(points_list, points)
            push!(equilibria_list, equilibria)
        end
    end
    codim_object(codim, points_list, equilibria_list)
end

function codim_map(ui, pi, codim::Tuple{Symbol, Symbol, Symbol};
        c1_lims = (-10.0, 10.0), c2_lims = (0.0, 10.0), c3_lims = (0.0, 10.0), resolution = 50, eq_res = 3,
    )
    uv = copy(ui);
    pv = copy(pi)
    points_list = Array{Tuple{Float64, Float64, Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], resolution)
    c2_range = LinRange(c2_lims[1], c2_lims[2], resolution)
    c3_range = LinRange(c3_lims[1], c3_lims[2], resolution)
    for (idx1, c1) in enumerate(c1_range)
        for (idx2, c2) in enumerate(c2_range)
            for (idx3, c3) in enumerate(c3_range)
                pv = copy(p0)
                uv = copy(u0)
                pv[(codim[1] |> p_find)...] = c1
                pv[(codim[2] |> p_find)...] = c2
                pv[(codim[3] |> p_find)...] = c3
                equilibria = find_equilibria(uv, pv; resolution = eq_res)
                points = (c1, c2, c3)
                #println(points |> typeof)
                push!(points_list, points)
                push!(equilibria_list, equilibria)
            end
        end
    end
    codim_object(codim, points_list, equilibria_list)
end

"""
A threshold is defined by the presence of a saddle node and a stable node. If only a stable node exists, spiking cannot occur, and if only a focus point exists, then spiking always occurs. 
"""
function extract_thresholds(c2::codim_object)
    if isa(c2.vars, Symbol)
        xs = []; threshs = [];
        for (idx, x) in enumerate(c2.points)
            eq = c2.equilibria[idx]
            if eq.stable != [] && eq.saddle != [] 
                push!(xs, x)
                push!(threshs, eq.saddle[1][1])
            end 
        end
        return xs, threshs          
    else
        xspan = unique(map(x -> x[1], c2.points))
        yspan = unique(map(x -> x[2], c2.points))
        nx = length(xspan); ny = length(yspan)
        thresh_map = zeros(nx, ny)
        for idx = 1:length(c2.points)
            eq = c2.equilibria[idx]
            if eq.stable != [] && eq.saddle == []
                thresh_map[idx] = Inf
            elseif eq.stable != [] && eq.saddle != [] 
                thresh_map[idx] = eq.saddle[1][1]
            elseif eq.unstable_focus != [] || eq.stable_focus != [] && eq.stable == []
                thresh_map[idx] = 0.0
            end
        end
        xspan, yspan, thresh_map
    end
end

"""
A baseline is defined by the presence of a stable node
"""
function extract_baselines(c2::codim_object)
    if isa(c2.vars, Symbol)
        xs = []; threshs = [];
        for (idx, x) in enumerate(c2.points)
            eq = c2.equilibria[idx]
            if eq.stable != []
                push!(xs, x)
                push!(threshs, eq.stable[1][1])
            elseif eq.stable_focus != [] 
                push!(xs, x)
                push!(threshs, eq.stable_focus[1][1])
            end
        end
        return xs, threshs
    else
        xspan = unique(map(x -> x[1], c2.points))
        yspan = unique(map(x -> x[2], c2.points))
        nx = length(xspan); ny = length(yspan)
        baseline_map = zeros(nx, ny)
        for idx = 1:length(c2.points)
            eq = c2.equilibria[idx]
            if eq.stable != [] 
                baseline_map[idx] = eq.stable[1][1]
            elseif eq.stable_focus != [] 
                baseline_map[idx] = eq.stable_focus[1][1]
            elseif eq.unstable_focus != [] && eq.stable == []
                baseline_map[idx] = Inf
            end
        end
        xspan, yspan, baseline_map
    end
end


