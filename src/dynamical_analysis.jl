"""
This function is a work of progress. Things I want it to include
- Using the symbols to declare the changing parameter or initial condition
- Changing multiple parameters at once
"""
function ensemble_func(prob::ODEProblem, i, repeat, idx, val_rng; run_func_on = :pars, verbose = false)
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

function ensemble_func(prob::SDEProblem, i, repeat, idx, val_rng; run_func_on = :pars, verbose = false)
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
        prob.u0[idx] .= val_rng[i]
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
    var_idxs = map(vr -> u_find(vr), vars)
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
        precision = 2, check_min = 1e-8
    )
    stable = Array{Array{Float64}}([])
    unstable = Array{Array{Float64}}([])
    saddle = Array{Array{Float64}}([])
    unstable_focus = Array{Array{Float64}}([])
    stable_focus = Array{Array{Float64}}([])
    storage = Array{Array{Float64}}([])

    var_idx = [(vars[1]|>u_find), (vars[2]|>u_find)]
    for (idx_x, x) in enumerate(LinRange(xlims[1], xlims[2], equilibrium_resolution)) 
        #Iterate through the x range
        for (idx_y, y) in enumerate(LinRange(ylims[1], ylims[2], equilibrium_resolution))
            #Iterate through the y range looking for stable points
            
            #We really may only need to check for stable points if the dU is low
            uI = prob.u0
            uI[var_idx] .= (x, y)
            df(ux) = prob.f(ux, prob.p, 0.0)
            res = nlsolve(df, uI)
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
                J_mat = ForwardDiff.jacobian(df, res.zero)[[var_idx...], [var_idx...]]
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
Make a 1D codim object
"""
function codim_map(prob, codim::Symbol;
        c1_lims = (-10.0, 10.0), codim_resolution = 50,
        saddle_continuation::Bool = true, 
        kwargs... #finding equilibria args
    )
    points_list = Array{Tuple{Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], codim_resolution)
    cont_toggle = false
    for (idx1, c1) in enumerate(c1_range)
        pv = prob.p
        pv[codim |> p_find] = c1
        prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
        equilibria = find_equilibria(prob_i; kwargs...)
        #in order to pass we want to make sure that there has at least been a saddle node first

        if cont_toggle == false && !isempty(equilibria.saddle)
            cont_toggle = true
            if saddle_continuation && length(points_list) > 1
                points_cont, eq_cont = eq_continuation(prob_i, (c1_range[idx1-1], c1_range[idx1]), codim; forward = false, kwargs...) #Step back
                splice!(points_list, idx1:idx1-1, points_cont)
                splice!(equilibria_list, idx1:idx1-1, eq_cont)
            end
        end
        
        if saddle_continuation && cont_toggle && isempty(equilibria.saddle)
            points_cont, eq_cont = eq_continuation(prob_i, (c1_range[idx1-1], c1_range[idx1]), codim) #Step back
            push!(points_list, points_cont...)
            push!(equilibria_list, eq_cont...)
            saddle_continuation = false #Use it only once
        end
        points = c1
        #println(points |> typeof)
        push!(points_list, (c1,))
        push!(equilibria_list, equilibria)
    end
    codim_object((codim,), points_list, equilibria_list)    
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
        #println("continuation: $(rng[1]) -> $(rng[2])")
        while I < In && ϵ > min_step && iter <= max_iters 
            iter += 1
            I += ϵ #Increment I slowly
            pv = prob.p
            pv[par |> p_find] = I #plug in the newly incremented equilibria
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
        #println("reverse continuation: $(rng[2]) -> $(rng[1])")
        while I > In && ϵ > min_step && iter < max_iters 
            iter += 1
            I -= ϵ #decrement I slowly
            pv = prob.p
            pv[par |> p_find] = I #plug in the newly incremented equilibria
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


function codim_map(prob::ODEProblem, codim::Tuple{Symbol, Symbol};
        c1_lims = (-10.0, 10.0), c2_lims = (0.0, 10.0), resolution = 50, eq_res = 3,
    )
    points_list = Array{Tuple{Float64, Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    c1_range = LinRange(c1_lims[1], c1_lims[2], resolution)
    c2_range = LinRange(c2_lims[1], c2_lims[2], resolution)
    for (idx1, c1) in enumerate(c1_range)
        #println(idx1)
        for (idx2, c2) in enumerate(c2_range)
            pv = prob.p
            pv[codim[1] |> p_find] = c1
            pv[codim[2] |> p_find] = c2
            prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
            equilibria = find_equilibria(prob_i; kwargs...)
            points = (c1, c2)
            #println(points |> typeof)
            push!(points_list, points)
            push!(equilibria_list, equilibria)
        end
    end
    codim_object(codim, points_list, equilibria_list)
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


