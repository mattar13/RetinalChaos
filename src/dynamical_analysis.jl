"""
This function is a work of progress. Things I want it to include
- Using the symbols to declare the changing parameter or initial condition
- Changing multiple parameters at once
"""

function ensemble_func(prob, i, repeat, idx::Int64, val_rng; run_func_on = :pars, verbose = false)
    if run_func_on == :pars
        if verbose
            println("Changing parameter $(prob.p[idx]) -> $(val_rng[i])")
        end
        prob.p[idx] .= val_rng[i]
        prob
    elseif run_func_on == :conds
        if verbose
            println("Changing condition $(prob.u0[idx]) -> $(val_rng[i])")
        end
        prob.u0 .= val_rng[i]
        prob
    end
end

function ensemble_func(prob, i, repeat, sym::Symbol, val_rng; run_func_on = :pars, verbose = false)
    idx = findall(x -> x==sym, tar_pars)
    if run_func_on == :pars
        if verbose
            println("Changing parameter $(prob.p[idx]) -> $(val_rng[i])")
        end
        prob.p[idx] .= val_rng[i]
        prob
    elseif run_func_on == :conds
        if verbose
            println("Changing condition $(prob.u0[idx]) -> $(val_rng[i])")
        end
        prob.u0 .= val_rng[i]
        prob
    end
end

function phase_plane(prob::ODEProblem; vars::Array{Symbol, 1} = [:v, :n], xlims = (-90.0, 10.0), ylims = (-0.10, 5.0), resolution = 100)
    var_idx = [(vars[1]|>u_find), (vars[2]|>u_find)]
    phase_plane = zeros(resolution, resolution, 2)
     for (idx_x, x) in enumerate(LinRange(xlims[1], xlims[2], resolution)), (idx_y, y) in enumerate(LinRange(ylims[1], ylims[2], resolution))
        uI = prob.u0
        uI[var_idx].= (x, y)
        du = prob.f(uI, prob.p, 0.0)
        phase_plane[idx_x, idx_y, :] = du[var_idx]
    end
    phase_plane
end

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

struct equilibria_object{T}
    stable::Array{Array{T}}
    unstable::Array{Array{T}}
    saddle::Array{Array{T}}
    unstable_focus::Array{Array{T}}
    stable_focus::Array{Array{T}}
end



#We can import a function that counts all of the equilibria events in the object
length(eq::equilibria_object) = length(eq.stable) + length(eq.unstable) + length(eq.saddle) + length(eq.unstable_focus) + length(eq.stable_focus)

#This function displays information about the equilibrium
function print(eq::equilibria_object; msg = "Existant Equilibrium", vars = [:v])
    println(msg)
    var_idxs = map(vr -> findall(x -> x==vr, tar_conds)[1], vars)
    eq.unstable != [] ? println("Unstable Equilibrium: $(eq.unstable[1][var_idxs])") : nothing
    eq.stable != [] ? println("Stable Equilibrium: $(eq.stable[1][var_idxs])") : nothing
    eq.saddle != [] ? println("Saddle Equilibrium: $(eq.saddle[1][var_idxs])") : nothing
    eq.unstable_focus != [] ? println("Unstable Focus Equilibrium: $(eq.unstable_focus[1][var_idxs])") : nothing
    eq.stable_focus != [] ? println("Stable Focus Equilibrium: $(eq.stable_focus[1][var_idxs])") : nothing
end

export print, length

#Conduct a stability analysis of the current
function find_equilibria(prob::ODEProblem;
        vars = [:v, :n], xlims = (-90.0, 10.0), ylims = (-1.0, 5.0), equilibrium_resolution = 10,
        precision = 2, check_min = 1e-5
    )
    stable = Array{Array{Float64}}([])
    unstable = Array{Array{Float64}}([])
    saddle = Array{Array{Float64}}([])
    unstable_focus = Array{Array{Float64}}([])
    stable_focus = Array{Array{Float64}}([])
    storage = Array{Array{Float64}}([])

    var_idx = [findall(x -> x==vars[1], tar_conds)[1], findall(x -> x==vars[2], tar_conds)[1]]
    for (idx_x, x) in enumerate(LinRange(xlims[1], xlims[2], equilibrium_resolution)) 
        #Iterate through the x range
        for (idx_y, y) in enumerate(LinRange(ylims[1], ylims[2], equilibrium_resolution))
            #Iterate through the y range looking for stable points
            
            #We really may only need to check for stable points if the dU is low
            uI = copy(prob.u0)
            uI[var_idx] .= (x, y)

            df(ux) = prob.f(similar(ux), ux, prob.p, 0.0) #dU, U, p, t            
            res = nlsolve(df, uI)
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
                J_mat = ForwardDiff.jacobian(df, res.zero)[[var_idx...], [var_idx...]]
                ev = eigvals(J_mat)

                if sign(real(ev[1])) != sign(real(ev[2]))
                    #If the sign of the eigenvalues are opposite, then it is a unstable saddle
                    println("Saddle")
                    push!(saddle, res.zero)
                else
                    if imag(ev[1]) != 0 #This detects whether or not the equilibria is a focus
                        if real(ev[1]) > 0.0
                            println("Unstable focus")
                            push!(unstable_focus, res.zero)
                        else
                            println("Unstable focus")
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

function eq_continuation(prob, rng::Tuple{T, T}, par::Symbol;
        forward = true, max_iters = 100, min_step = 1.0e-15, 
        kwargs...
    ) where T <: Real

    points_list = Array{Tuple{Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    bifurcation_point = 0.0
    bifurcation_eq = nothing
    #Adaptive continuation
    if forward
        I = rng[1] #Step back
        In = rng[2] #we start here 
        ϵ = abs(I - In)/2 #Begin at the halfway point between the two points
        iter = 0
        println("continuation: $(rng[1]) -> $(rng[2])")
        while I < In && ϵ > min_step && iter <= max_iters 
            iter += 1
            I += ϵ #Increment I slowly
            pv = prob.p
            pv[findall(x -> x==par, tar_pars)[1]] = I #plug in the newly incremented equilibria
            prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
            equilibria = find_equilibria(prob_i; kwargs...)

            #println(points |> typeof)
            #If the number of saddle equilibria drops to 0, then return to the previous
            if length(equilibria) == 2 
                #This is a flaw of the find equilibria algorithim, move away from this point
                I -= ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            elseif isempty(equilibria.saddle) #The saddle node is terminated past here
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
                I -= ϵ #walk back
                ϵ /= 2 #Divide epsilon in half

            else #None of these contiditons were met alter the bifurcation eq
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
            end

        end
    else #forward == false
        I = rng[2] #Step back
        In = rng[1] #we start here 
        ϵ = abs(I - In)/2 #Begin at the halfway point between the two points
        iter = 0
        println("reverse continuation: $(rng[2]) -> $(rng[1])")
        while I > In && ϵ > min_step && iter < max_iters 
            iter += 1
            I -= ϵ #decrement I slowly
            pv = prob.p
            pv[findall(x -> x==par, tar_pars)[1]] = I #plug in the newly incremented equilibria
            prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
            equilibria = find_equilibria(prob_i; kwargs...)
            #println(points |> typeof)
            #we will add each point and new equilibria to the solution
            #If the number of saddle equilibria drops to 0, then return to the previous
              if length(equilibria) == 2 
                #This is a flaw of the find equilibria algorithim, move away from this point
                I += ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            elseif isempty(equilibria.saddle) #The saddle node is terminated past here
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
                I += ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            else #None of these contiditons were met alter the bifurcation eq
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
            end
        end
    end
    #we want to sort the points before adding them
    points_vals = map(point -> point[1], points_list)
    sort_idxs = sortperm(points_vals)
    #println(sort_idxs)
    return points_list[sort_idxs], equilibria_list[sort_idxs]
end

"""
Make a 1D codim object
"""
function codim_map(prob, codim::Symbol, c1_lims::Tuple{T, T};
        codim_resolution = 50, saddle_continuation::Bool = true, 
        kwargs... #finding equilibria args
    ) where T <: Real
    points_list = Array{Tuple{T}}([])
    equilibria_list = Array{equilibria_object{T}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], codim_resolution)
    cont_toggle = false
    n_equilibria = -1
    for (idx1, c1) in enumerate(c1_range)
        println(idx1)
        pv = prob.p
        pv[findall(x -> x==codim, tar_pars)[1]] = c1
        #println(pv[findall(x -> x==codim, tar_pars)[1]])
        #println(pv)
        #println(tar_pars[findall(x -> x==codim, tar_pars)[1]])
        prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
        equilibria = find_equilibria(prob_i; kwargs...)
        #in order to pass we want to make sure that there has at least been a saddle node first
        
        if n_equilibria == -1 #This is the starting point
            n_equilibria = length(equilibria)
        elseif length(equilibria) > n_equilibria
            if saddle_continuation
                points_cont, eq_cont = eq_continuation(prob_i, (c1_range[idx1-1], c1_range[idx1]), codim; forward = false, kwargs...) #Step back
                splice!(points_list, idx1:idx1-1, points_cont)
                splice!(equilibria_list, idx1:idx1-1, eq_cont)
            end
            n_equilibria = length(equilibria)
        elseif length(equilibria) < n_equilibria
            if saddle_continuation
                points_cont, eq_cont = eq_continuation(prob_i, (c1_range[idx1-1], c1_range[idx1]), codim) #Step back
                push!(points_list, points_cont...)
                push!(equilibria_list, eq_cont...)
            end
            #saddle_continuation = false #Use it only once
            n_equilibria = length(equilibria)
        elseif length(equilibria) == n_equilibria
            n_equilibria = length(equilibria)
        end
        print
        points = c1
        push!(points_list, (c1,))
        push!(equilibria_list, equilibria)
    end
    codim_object((codim,), points_list, equilibria_list)    
end

#This is for a 1 dimensional array
"""
This function finds a bifurcation which is the elimination of a saddle node at a junction
"""
function find_bifurcation(c1_map::codim_object{1, T}) where T <: Real
    c_vals = map(x -> x[1], c1_map.points)
    bif_value = T[]
    bif_eq = []

    #println(n_eqs)
    saddle_rng = findall(map(x -> length(x.saddle) >= 1, c1_map.equilibria))
    if saddle_rng[1] != 1
        push!(bif_value, c1_map.points[saddle_rng[1]][1])
        push!(bif_eq, c1_map.equilibria[saddle_rng[1]])
    end
    
    if saddle_rng[end] != length(c_vals)
        push!(bif_value, c1_map.points[saddle_rng[end]][1])
        push!(bif_eq, c1_map.equilibria[saddle_rng[end]])
    end
    #Method 2 for finding bifurcations
    #n_eqs = map(x -> length(x), c1_map.equilibria)
    #change = n_eqs[1:end-1] .- n_eqs[2:end]
    #pre_saddle = findall(change.==-2).+1 #detect a decrease in equilibrium by 2
    #post_saddle = findall(change.==2) #detect a decrease in equilibrium by 2

    #Push bifurcations for reverse continuation
    #if !isempty(pre_saddle)
    #    push!(bif_value, c1_map.points[pre_saddle][1][1])
    #    push!(bif_eq, c1_map.equilibria[pre_saddle][1])
    #end
    #Push bifurcations for continuation
    #if !isempty(post_saddle)
    #    push!(bif_value, c1_map.points[post_saddle][1][1])
    #    push!(bif_eq, c1_map.equilibria[post_saddle][1])
    #end
    return bif_value, bif_eq
end

function codim_map(prob::ODEProblem, codim::Tuple{Symbol, Symbol}, c1_lims::Tuple{T, T}, c2_lims::Tuple{T, T};  
        codim_resolution = 10,
        kwargs...
    ) where T <: Real
    points_list = Array{Tuple{T, T}}([])
    equilibria_list = Array{equilibria_object{T}}([])
    c2_range = LinRange(c2_lims[1], c2_lims[2], codim_resolution)
    #for (idx1, c1) in enumerate(c1_range) #We will do the equilibria continuation on this parameter
    for (idx2, c2) in enumerate(c2_range)
        pv = prob.p
        #pv[codim[1] |> p_find] = c1
        pv[codim[2] |> p_find] = c2
        prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
        #equilibria = find_equilibria(prob_i; kwargs...)
        #Conduct an inner loop equilibria with continution focused on param c1
        c_inner = codim_map(prob_i, codim[1], c1_lims, codim_resolution = codim_resolution)
        points = map(x -> (x[1],c2), c_inner.points)
        push!(points_list, points...)
        push!(equilibria_list, c_inner.equilibria...)
    end
    codim_object(codim, points_list, equilibria_list)
end

#This is for a 2 dimensional array
function find_bifurcation(c2_map::codim_object{2, T}; idx::Int64 = 1) where T <: Real
    c_vals_c1 = map(x -> x[1], c2_map.points)
    c_vals_c2 = map(x -> x[2], c2_map.points)
    println(unique(c_vals_c1))
    #saddle_rng = findall(map(x -> length(x.saddle) >= 1, c1_map.equilibria))
    #bif_value = c2_map.points[saddle_rng[end]][1]
    #bif_eq = c2_map.equilibria[saddle_rng[end]].saddle
    #return bif_value, bif_eq
end

function codim_map(prob::ODEProblem, codim::Tuple{Symbol, Symbol, Symbol};
        c1_lims = (-10.0, 10.0), c2_lims = (0.0, 10.0), c3_lims = (0.0, 10.0), resolution = 50, eq_res = 3,
    )
    points_list = Array{Tuple{Float64, Float64, Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], resolution)
    c2_range = LinRange(c2_lims[1], c2_lims[2], resolution)
    c3_range = LinRange(c3_lims[1], c3_lims[2], resolution)
    for (idx1, c1) in enumerate(c1_range)
        for (idx2, c2) in enumerate(c2_range)
            for (idx3, c3) in enumerate(c3_range)
                pv = prob.p
                pv[(codim[1] |> p_find)...] = c1
                pv[(codim[2] |> p_find)...] = c2
                pv[(codim[3] |> p_find)...] = c3
                prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
                equilibria = find_equilibria(prob_i; resolution = eq_res)
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


