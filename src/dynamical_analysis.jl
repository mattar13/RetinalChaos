norminf(x) = norm(x, Inf)

"""
ensemble_func()  sets up a ensemble problem
"""

function ensemble_func(prob, i, repeat, idx::Int64, val_rng; verbose = true)
    #println("Running this version")
    if verbose
        println("Changing parameter $(prob.p[idx]) -> $(val_rng[i])")
    end
    prob.p[idx] = val_rng[i]
    prob
end

ensemble_func(prob, i, repeat, sym::Symbol, val_rng; verbose=true) = ensemble_func(prob, i , repeat, indexof(sym), val_rng; verbose = verbose)

"""
monte_func() sets up a monte carlo problem
"""
function monte_func(prob, i, repeat; pars = :all, dists = [Normal()])
    if pars == :all
        new_pars = map(x -> typeof(x) != Float64 ? rand(x) : x , dists);
        if new_pars[(:σ |> p_find)[1]] > 0.0
            stochastic = true
        end
        #new_conds = map(x -> typeof(x) != Float64 ? rand(x) : x,  dists);
        if stochastic
            return SDEProblem(BurstModel, noise, prob.u0, prob.tspan, new_pars)
        else
            return ODEProblem(BurstModel, prob.u0, prob.tspan, new_pars)
        end
    else
        stochastic = false
        new_pars = prob.p
        new_conds = prob.u0
        if prob.p[(:σ |> p_find)[1]] > 0.0
            stochastic = true
        end
        for (idx, var) in enumerate(pars)
            conds = var |> u_find
            pars = var |> p_find
            if length(conds) > 0
                new_conds[conds[1]] = rand(dists[idx])
            elseif length(pars) > 0
                new_pars[pars[1]] = rand(dists[idx])
            end
        end

        if stochastic
            return SDEProblem(BurstModel, noise, new_conds, prob.tspan, new_pars)
        else
            return ODEProblem(BurstModel, new_conds, prob.tspan, new_pars)
        end
    end
end


"""
Current_Clamp experiment

In order to conduct one of these experiments we can set up an epoch table. 
"""
function IC_callback(step_begin, duration, level; duration_delta = 0, level_delta = 0, give_warning = false)
    if duration_delta == 0 && level_delta == 0 #We are only doing one duration vs a tar_conds
        condition_fn1(u, t, integrator) = step_begin < t < step_begin + duration
        affect_fn1!(integrator) = integrator.p[:I_app |> p_find] = level
        cb_fn1 = DiscreteCallback(condition_fn1, affect_fn1!)
        return cb_fn1
    elseif duration_delta != 0 && level_delta == 0
        #These callbacks will be called in an ensemble solution
        if give_warning
            @warn begin
                "These solutions will instead return an ensemble solution"
            end
        end
        condition_fn2(u, t, integrator, i) = step_begin < t < step_begin + (duration + duration_delta*i)
        affect_fn2!(integrator) = integrator.p[:I_app |> p_find] = level
        cb_fn2(i) = DiscreteCallback((u, t, integrator) -> condition_fn2(u, t, integrator,  i), affect_fn2!)
        return cb_fn2
    elseif duration_delta != 0 && level_delta != 0
        if give_warning
            @warn begin
                "These solutions will only work with an ensemble solution"
            end
        end
        condition_fn3(u, t, integrator) = step_begin < t < step_begin + duration
        affect_fn3!(integrator, i) = integrator.p[:I_app |> p_find] = level + level_delta*i
        cb_fn3(i) = DiscreteCallback(condition_fn3, integrator -> affect_fn3!(integrator, i))
        return cb_fn3
    else
        if give_warning
            @warn begin
                "These solutions will only work with an ensemble solution"
            end
        end
        condition_fn4(u, t, integrator, i) = step_begin < t < step_begin + (duration + duration_delta*i)
        affect_fn4!(integrator, i) = integrator.p[:I_app |> p_find] = level + (level_delta*i)
        cb_fn4(i) = DiscreteCallback(
            (u, t, integrator) -> condition_fn4(u, t, integrator,  i), 
            integrator -> affect_fn4!(integrator, i)
            )
        return cb_fn4
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
    stable_focus::Array{Array{T}}
    unstable_focus::Array{Array{T}}
end

import Base.getindex
function getindex(eq::equilibria_object, sym::Symbol)
    if !isnothing(sym |> u_find)
        var = sym |> u_find
        a = [
            !isempty(eq.stable) ? eq.stable[1][var] : NaN,
            !isempty(eq.unstable) ? eq.unstable[1][var] : NaN,
            !isempty(eq.saddle) ? eq.saddle[1][var] : NaN,
            !isempty(eq.stable_focus) ? eq.stable_focus[1][var] : NaN,
            !isempty(eq.unstable_focus) ? eq.unstable_focus[1][var] : NaN
        ]
        return a
    elseif sym == :stable
        return eq.stable
    elseif sym == :unstable
        return eq.unstable
    elseif sym == :saddle
        return eq.saddle
    elseif sym == :stable_focus
        return eq.stable_focus
    elseif sym == :unstable_focus
        return eq.unstable_focus
    end
end

getindex(eq::equilibria_object, syms...) = map(sym -> eq[sym], syms)

function getindex(eq::equilibria_object, idx::Int64)  
    if idx == 1
        return eq.stable
    elseif idx == 2
        return eq.unstable
    elseif idx == 3
        return eq.saddle
    elseif idx == 4
        return eq.stable_focus
    elseif idx == 5
        return eq.unstable_focus
    end
end

#We can import a function that counts all of the equilibria events in the object
length(eq::equilibria_object) = length(eq.stable) + length(eq.unstable) + length(eq.saddle) + length(eq.unstable_focus) + length(eq.stable_focus)

#This function displays information about the equilibrium
function print(eq::equilibria_object; msg = "Existant Equilibrium", vars = [:v])
    println(msg)
    #var_idxs = map(vr -> findall(x -> x==vr, tar_conds)[1], vars)
    var_idxs = map(vr -> vr |> u_find, vars)
    eq.unstable != [] ? println("Unstable Equilibrium: $(eq.unstable[1][var_idxs])") : nothing
    eq.stable != [] ? println("Stable Equilibrium: $(eq.stable[1][var_idxs])") : nothing
    eq.saddle != [] ? println("Saddle Equilibrium: $(eq.saddle[1][var_idxs])") : nothing
    eq.unstable_focus != [] ? println("Unstable Focus Equilibrium: $(eq.unstable_focus[1][var_idxs])") : nothing
    eq.stable_focus != [] ? println("Stable Focus Equilibrium: $(eq.stable_focus[1][var_idxs])") : nothing
end

export print, length

#Conduct a stability analysis of the current
#We have to ensure the points we pick are realistic
function find_equilibria(prob::ODEProblem;
        vars = [:v, :n], 
        xlims = (-100.0, 10.0), ylims = (0.0, 1.0), #These need to be in realistic ranges. 
        equilibrium_resolution = 10,
        precision = 2, check_min = 1e-5, verbose = false
    )
    stable = Array{Array{Float64}}([])
    unstable = Array{Array{Float64}}([])
    saddle = Array{Array{Float64}}([])
    unstable_focus = Array{Array{Float64}}([])
    stable_focus = Array{Array{Float64}}([])
    storage = Array{Array{Float64}}([])

    #var_idx = [findall(x -> x==vars[1], tar_conds)[1], findall(x -> x==vars[2], tar_conds)[1]]
    var_idx = [vars[1] |> u_find, vars[2] |> u_find]
    for (idx_x, x) in enumerate(LinRange(xlims[1], xlims[2], equilibrium_resolution)) 
        #Iterate through the x range
        for (idx_y, y) in enumerate(LinRange(ylims[1], ylims[2], equilibrium_resolution))
            #Iterate through the y range looking for stable points
            if verbose
                println("Checking parameter")
                println("$(vars[1]) = $x")
                println("$(vars[2]) = $y")
            end
            uI = copy(prob.u0)
            uI[var_idx] .= (x, y)
        
            if verbose
                println("Condition to check = $uI")
            end
        
            df(ux) = prob.f(similar(ux), ux, prob.p, 0.0) #dU, U, p, t            
            res = nlsolve(df, uI)
            #If any of the results are out of bounds then we have to cut it off
            if verbose
                print("Results of the nlsolve for uI: ")
                #println(res)
                print("f = 0: ")
                println(res.zero)
            end
        
            #We want to remove all cases where the result is out of bounds
            var1_inbounds = xlims[1] < res.zero[var_idx[1]] <= xlims[2]
            var2_inbounds = ylims[1] < res.zero[var_idx[2]] <= ylims[2]

            if var1_inbounds && var2_inbounds #We want to make sure the equilibria is not out of bounds
                equilibria = map(x -> round(x, digits=precision), res.zero)
                check = df(res.zero)
                if verbose
                    print("Checking the solution and it's proximity to zero: dUᵢ = ")
                    println(check)
                end
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
                        #println("Saddle")
                        push!(saddle, res.zero)
                    else
                        if imag(ev[1]) != 0 #This detects whether or not the equilibria is a focus
                            if real(ev[1]) >= 0.0
                                #println("Unstable focus")
                                push!(unstable_focus, res.zero)
                            else
                                #println("Unstable focus")
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
            else
                if !var1_inbounds && verbose
                    println("Variable $(vars[1]) is out of bounds at $(res.zero[var_idx[1]]) will be skipped")
                end

                if !var2_inbounds && verbose
                    println("Variable $(vars[2]) is out of bounds at $(res.zero[var_idx[2]]) will be skipped")
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

function getindex(c1::codim_object{1, T}, sym::Symbol) where T <: Real
    
    if !isnothing(sym |> u_find) #This will occur if the symbol is a ic symbol
        var_idx = sym |> u_find
        return hcat(
            map(eq -> length(eq.saddle) > 0 ? eq.saddle[1][var_idx] : NaN, c1.equilibria),
            map(eq -> length(eq.stable) > 0 ? eq.stable[1][var_idx] : NaN, c1.equilibria),
            map(eq -> length(eq.unstable) > 0 ? eq.unstable[1][var_idx] : NaN, c1.equilibria),
            map(eq -> length(eq.unstable_focus) > 0 ? eq.unstable_focus[1][var_idx] : NaN, c1.equilibria),
            map(eq -> length(eq.stable_focus) > 0 ? eq.stable_focus[1][var_idx] : NaN, c1.equilibria)
        )
        
    elseif sym == :points
        points = map(x -> x[1], c1.points);
    end
end

getindex(c1::codim_object{1, T}, syms...) where T <: Real = map(sym -> c1[sym], syms)

function eq_continuation(prob, rng::Tuple{T, T}, par::Symbol;
        forward = true, max_iters = 100, min_step = 1.0e-10, 
        record_noisy_points = true,
        kwargs...
    ) where T <: Real

    points_list = Array{Tuple{Float64}}([])
    equilibria_list = Array{equilibria_object{Float64}}([])
    #Adaptive continuation
    if forward
        I = rng[1] #Step back
        In = rng[2] #we start here 
        ϵ = abs(I - In) / 2 #Begin at the halfway point between the two points
        iter = 0
        #println("$par continuation: $(rng[1]) -> $(rng[2])")
        while I < In && ϵ > min_step && iter <= max_iters
            iter += 1
            I += ϵ #Increment I slowly
            pv = prob.p
            pv[par|>p_find] = I #plug in the newly incremented equilibria
            prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
            equilibria = find_equilibria(prob_i; kwargs...)
            #println("$par = $I results in $(length(equilibria)) equilibria")
            #println("Epsilon at $ϵ")
            #println("Iterations at $(iter)")
            #If the number of saddle equilibria drops to 0, then return to the previous
            if length(equilibria) == 2
                #This is a flaw of the find equilibria algorithim, move away from this point
                #println("$par has found a noisy equilibria pair at $I (usually indicating near annhilation)")
                if record_noisy_points
                    push!(points_list, (I,))
                    push!(equilibria_list, equilibria)
                end
                I -= ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            elseif isempty(equilibria.saddle) #The saddle node is terminated past here
                #println("$par has found a saddle node at $I")
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
                I -= ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            
            else #None of these contiditons were met alter the bifurcation eq
                #println("$par has not yet found a equilibria at $I")
                #I += ϵ #Adding in another leap forward may speed things up a bit
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
            end
    
        end
    else #forward == false
        I = rng[2] #Step back
        In = rng[1] #we start here 
        ϵ = abs(I - In) / 2 #Begin at the halfway point between the two points
        iter = 0
        #println("$par reverse continuation: $(rng[1]) <- $(rng[2])")
        while I > In && ϵ > min_step && iter < max_iters
            iter += 1
            I -= ϵ #decrement I slowly
            pv = prob.p
            pv[par|>p_find] = I #plug in the newly incremented equilibria
            prob_i = ODEProblem(prob.f, prob.u0, prob.tspan, pv)
            equilibria = find_equilibria(prob_i; kwargs...)
            #println("$par = $I results in $(length(equilibria)) equilibria")
            #println("Epsilon at $ϵ")
            #println("Iterations at $(iter)")
            #we will add each point and new equilibria to the solution
            #If the number of saddle equilibria drops to 0, then return to the previous
            if length(equilibria) == 2
                #println("$par has found a noisy equilibria pair at $I (usually indicating near annhilation)")
                if record_noisy_points
                    push!(points_list, (I,))
                    push!(equilibria_list, equilibria)
                end
                #This is a flaw of the find equilibria algorithim, move away from this point
                I += ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            elseif isempty(equilibria.saddle) #The saddle node is terminated past here
                #println("$par reverse has found a saddle node at $I")
                push!(points_list, (I,))
                push!(equilibria_list, equilibria)
                I += ϵ #walk back
                ϵ /= 2 #Divide epsilon in half
            else #None of these contiditons were met alter the bifurcation eq
                #println("$par reverse has not yet found a equilibria at $I")
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
    n_equilibria = -1
    for (idx1, c1) in enumerate(c1_range)
        #println("Setting variable $codim $c1")
        pv = prob.p
        #pv[findall(x -> x==codim, tar_pars)[1]] = c1
        pv[codim |> p_find] = c1
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
        println("$(n_equilibria) found at variable $codim $c1")
        #points = c1
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

function segment_extracter(c1::codim_object{1, T}, vars = :v) where T <: Real
    points = map(x -> x[1], c1.points);
    equilibria = zeros()
    saddle_p = map(eq -> length(eq.saddle) > 0 ? eq.saddle[1][var_idx] : NaN, c1.equilibria);
    stable_p = map(eq -> length(eq.stable) > 0 ? eq.stable[1][var_idx] : NaN, c1.equilibria);
    unstable_p = map(eq -> length(eq.unstable) > 0 ? eq.unstable[1][var_idx] : NaN, c1.equilibria);
    unstable_focus_p = map(eq -> length(eq.unstable_focus) > 0 ? eq.unstable_focus[1][var_idx] : NaN, c1.equilibria);
    stable_focus_p = map(eq -> length(eq.stable_focus) > 0 ? eq.stable_focus[1][var_idx] : NaN, c1.equilibria);
    
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


