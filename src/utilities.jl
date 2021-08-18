"""
ensemble_func()  sets up a ensemble problem
"""

function ensemble_func(prob, i, repeat, idx::Int64, val_rng; run_func_on = :pars, verbose = false)
    if run_func_on == :pars
        if verbose
            println("Changing parameter $(prob.p[idx]) -> $(val_rng[i])")
        end
        prob.p[idx] = val_rng[i]
        prob
    elseif run_func_on == :conds
        if verbose
            println("Changing condition $(prob.u0[idx]) -> $(val_rng[i])")
        end
        prob.u0[idx] = val_rng[i]
        prob
    end
end

function ensemble_func(prob, i, repeat, sym::Symbol, val_rng; verbose = false)
    idx_cond = sym |> u_find
    idx_par = sym |> p_find

    if !isnothing(idx_par) 
        if verbose
            println("Changing parameter $(prob.p[idx]) -> $(val_rng[i])")
        end
        prob.p[idx_par] = val_rng[i]
    end    
    
    if !isnothing(idx_cond)
        if verbose
            println("Changing condition $(prob.u0[idx]) -> $(val_rng[i])")
        end
        prob.u0[idx_cond] = val_rng[i]
    end
    
    prob
end

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
        condition(u, t, integrator) = step_begin < t < step_begin + duration
        affect!(integrator) = integrator.p[:I_app |> p_find] = level
        cb = DiscreteCallback(condition, affect!)
        return cb
    elseif duration_delta != 0 && level_delta == 0
        #These callbacks will be called in an ensemble solution
        if give_warning
            @warn begin
                "These solutions will instead return an ensemble solution"
            end
        end
        condition(u, t, integrator, i) = step_begin < t < step_begin + (duration + duration_delta*i)
        affect!(integrator) = integrator.p[:I_app |> p_find] = level
        cb(i) = DiscreteCallback((u, t, integrator) -> condition(u, t, integrator,  i), affect!)
        return cb
    elseif duration_delta != 0 && level_delta != 0
        if give_warning
            @warn begin
                "These solutions will only work with an ensemble solution"
            end
        end
        condition(u, t, integrator) = step_begin < t < step_begin + duration
        affect!(integrator, i) = integrator.p[:I_app |> p_find] = level + level_delta*i
        cb(i) = DiscreteCallback(condition, integrator -> affect!(integrator, i))
        return cb
    else
        if give_warning
            @warn begin
                "These solutions will only work with an ensemble solution"
            end
        end
        condition(u, t, integrator, i) = step_begin < t < step_begin + (duration + duration_delta*i)
        affect!(integrator, i) = integrator.p[:I_app |> p_find] = level + (level_delta*i)
        cb(i) = DiscreteCallback(
            (u, t, integrator) -> condition(u, t, integrator,  i), 
            integrator -> affect!(integrator, i)
            )
        return cb
    end    
end
#Reading and writing JSON files
"""
This file writes either a named tuple or a dictionary into a JSON file
"""
function write_JSON(data::T, filename::String) where T
    string_data = JSON2.write(data)
    open(filename, "w") do f
        write(f, string_data)
    end
end

"""
This function reads JSON files into whatever type is_type is.
"""
function read_JSON(::Type{T}, name_file::String) where T <: Dict
    nt = nothing
    open(name_file, "r") do f
        nt = JSON2.read(f, T)
    end
    nt
end

read_JSON(name_file::String) = read_JSON(Dict{Symbol, Float64}, name_file)
"""
Extract a parameter, condition dictionary
USAGE:
#p_dict -> [The parameter dictionary read directly from JSON]
> p0 = extract_dict(p_dict, BurstModel.params)

If your initial conditions are for a network, you can add the extra argument for dims

> In[1]: dims = (64,64);
> In[2]: p0 = extract_dict(p_dict, BurstModel.params, dims)
> Out[2]: 64x64xN AbstractArray:
"""
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}) = map(x -> Float64(dict_item[x]), pars)
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}, dims::Tuple) = cat(map(x -> fill(Float64(dict_item[x]), dims), pars)..., dims = length(dims)+1)

#get the index of the conditions in the list
function u_find(cond::Symbol; list_u::Array{Symbol, 1} = tar_conds, safe_skip::Bool = false) 
    try
        return findall(c -> c == cond, list_u)[1]
    catch
        if safe_skip
            println("Value does not exist in the range")
        end
        return nothing
    end
end

#get the index of the parameter in the list
function p_find(par::Symbol; list_p::Array{Symbol, 1}  = tar_pars, safe_skip::Bool = false) 
    try
        findall(p -> p == par, list_p)[1]
    catch
        if safe_skip
            println("Value does not exist in the range")
        end
        return nothing
    end
end

function extract_dict(dict_item::Dict{Symbol, T}) where T <: Real
    #To extract the dictionary first we want to take a look at the first param to see if it is a parameter or a condition
    if dict_item.keys[1] ∈ tar_pars
        #This is a parameter and we should extract it
        pars = T[]
        for par in tar_pars
            push!(pars, dict_item[par])
        end
        return pars
    else
        conds = T[]
        for cond in tar_conds
            push!(conds, dict_item[cond])
        end
        return conds
    end
end

function extract_dict(dict_item::Dict{Symbol, T}, dims...) where T <: Real
    #Extract the dict_item
    vals = extract_dict(dict_item)
    dims = Int64.(dims) #Put this in to convert the dimensions  
    val_map = zeros(dims..., length(vals))
    for i in 1:length(vals)
        val_map[:,:,i] .= vals[i]
    end
    val_map
end

#This function converts the solution into a CPU based solution
function convert_to_cpu(model::RODESolution)
    prob_remade = remake(model.prob, u0 = model.prob.u0 |> Array)
    #println(prob_remade |> typeof)
    #lets go into remaking the problem
    alg = model.alg
    t = model.t
    u = map(u -> u |> Array, model.u) #This pa
    W = model.W |> Array
    return RetinalChaos.DiffEqBase.build_solution(prob_remade, alg, t, u, W=W)
end

"""
This function saves the solution. 
In order to save the solution correctly, ensure to run: 
    julia> sol = sol |> convert_to_cpu
    julia> save_solution(sol, save_path)

"""
function save_solution(sol, save_path::String; name = "sol", partitions = 1, mode = :bson)
    
    if mode == :bson
        if partitions == 1
            file_contents = Dict(:t => ["$(name)_t.bson"], :u => ["$(name)_u.bson"])
            #write the details to a 
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            bson("$(save_path)\\$(name)_t.bson", Dict(:t => sol.t))
            bson("$(save_path)\\$(name)_u.bson", Dict(:u => sol.u))
        else
            file_contents = Dict(:t => ["$(name)1_t.bson"], :u => ["$(name)1_u.bson"])
            for i in 2:partitions
                push!(file_contents[:t], "$(name)$(i)_t.bson")
                push!(file_contents[:u], "$(name)$(i)_u.bson")
            end
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            partition_idxs = round.(Int64, LinRange(1, length(sol.t), partitions+1))
            for idx in 1:length(partition_idxs)-1
                if idx == 1
                    bson("$(save_path)\\$(name)$(idx)_t.bson", Dict(:sol_t => sol.t[1:partition_idxs[idx+1]]))
                    bson("$(save_path)\\$(name)$(idx)_u.bson", Dict(:sol_u => sol.u[1:partition_idxs[idx+1]]))
                else
                    bson("$(save_path)\\$(name)$(idx)_t.bson", Dict(:sol_t => sol.t[partition_idxs[idx]+1:partition_idxs[idx+1]]))
                    bson("$(save_path)\\$(name)$(idx)_u.bson", Dict(:sol_u => sol.u[partition_idxs[idx]+1:partition_idxs[idx+1]]))
                end
            end
        end
    else
        println("TODO implement JLD2")
    end
end

function load_solution(load_path)
    file_contents = read_JSON(Dict{Symbol, Vector{String}}, "$(load_path)\\file_contents.json")
    n_partitions = length(file_contents[:u])
    if n_partitions == 1 #there are no partitions in the file
        sol_t = BSON.load("$(load_path)\\$(file_contents[:t][1])")[:sol_t]
        sol_u = BSON.load("$(load_path)\\$(file_contents[:u][1])")[:sol_u]
    else
        sol_t_arr  = Float64[]
        sol_u_arr = Vector{Vector{Float64}}()
        for i in 1:n_partitions
            sol_t = BSON.load("$(load_path)\\$(file_contents[:t][i])")[:sol_t]
            sol_u = BSON.load("$(load_path)\\$(file_contents[:u][i])")[:sol_u]
            push!(sol_t_arr, sol_t...)
            push!(sol_u_arr, sol_u...)
            println(size(sol_u_arr))
        end
    end

    #catch
    #    println("Loading method 2")
    #    sol_array_pt1 = BSON.load("$(load_path)\\sol_array_pt1.bson")[:sol_u]
    #    sol_array_pt2 = BSON.load("$(load_path)\\sol_array_pt2.bson")[:sol_u]
    #    vcat(sol_array_pt1, sol_array_pt2)
    #end
    #p_dict = read_JSON(Dict{Symbol, Float32}, "$(load_path)\\params.json")
    #p0 = p_dict |> extract_dict
    #net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ])
    #sol_prob = SDEProblem(net, noise, sol_u[1], (0f0 , p_dict[:t_warm]), p0)
    #SciMLBase.build_solution(sol_prob, SOSRI(), sol_t, sol_u) #we can use this to build a solution without GPU
end

"""
This function runs the model using the indicated parameters

"""
function run_model(file_root::String, p_dict::Dict{Symbol, T}, u_dict::Dict{Symbol, T}; 
            gpu::Bool = true, version = :gACh,
            abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
            save_sol = true, save_partitions = 1,
            animate_solution = true, animate_dt = 60.0, 
            model_file_type = :bson, 
            iterations = 1 #this option sets the model up into sections so that we can break up the saving of the solutions
        ) where T <: Real
    
    #First write the parameters to the root
    write_JSON(p_dict, "$(file_root)\\params.json") #write the parameters to save for later
    write_JSON(u_dict, "$(file_root)\\iconds.json") #write the 
    
    #Load the network interface
    if gpu
        u0 = extract_dict(u_dict, p_dict[:nx], p_dict[:ny]) |> cu
        p0 = p_dict |> extract_dict
    else
        u0 = extract_dict(u_dict, p_dict[:nx], p_dict[:ny])
        p0 = p_dict |> extract_dict
    end

    #Load the model and ODE interface
    net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ], version = version, gpu = gpu)
    NetProb = SDEProblem(net, noise, u0, (0f0 , p_dict[:t_warm]), p0)

    print("[$(now())]: Warming up the solution... ")
    @time sol = solve(NetProb, SOSRI(), 
        abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
        progress = true, progress_steps = 1, 
        save_everystep = false
    )

    warmup = sol[end] #Keep this one as the original 
    #Save the warmed up solution
    if model_file_type == :bson #In this case we want to make a backup of the file
        BSON.@save "$(file_root)\\conds.bson" warmup #as a BSON file
    elseif model_file_type == :jld2  #In this case we want to make a backup of the file
        JLD2.@save "$(file_root)\\conds.jld2" warmup #This is only here to try to save the older files
    end

    #Now we want to run the actual simulation
    if iterations == 1
        print("[$(now())]: Running the model... ")
        NetProb = SDEProblem(net, noise, warmup, (0f0 , p_dict[:t_run]), p0)
        #Run the solution to fruition
        @time sol = solve(NetProb, SOSRI(), 
            abstol = abstol, reltol = reltol, maxiters = maxiters,
            save_idxs = [1:(Int64(p_dict[:nx]*p_dict[:ny]))...], 
            progress = true, progress_steps = 1
        )
        sol = convert_to_cpu(sol) #Before saving we need to bump the file to CPU
        if save_sol
            print("[$(now())]: Saving the simulation...")
            save_solution(sol, file_root; mode = model_file_type, partitions = save_partitions) 
            println("Completed")
        end
    else
        #for i in iterations
    end
    
    if animate_solution
        println("[$(now())]: Animating simulation...")
        anim = @animate for t = 1.0:animate_dt:sol.t[end]
            frame_i = reshape(sol(t) |> Array, (p_dict[:nx]|>Int64, p_dict[:ny]|>Int64))
            heatmap(frame_i, ratio = :equal, grid = false,
                    xaxis = "", yaxis = "", xlims = (0, p_dict[:nx]), ylims = (0, p_dict[:ny]),
                    c = :curl, clims = (-70.0, 0.0),
            )
        end
        gif(anim, "$(file_root)\\animation.gif", fps = 1000.0/animate_dt)
    end

    #We want to return the solution
    return sol
end 