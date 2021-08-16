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

#Reading and writing JSON files
"""
This file writes either a named tuple or a dictionary into a JSON file
"""
function write_JSON(data, name_file)
    string_data = JSON2.write(data)
    open(name_file, "w") do f
        write(f, string_data)
    end
end

"""
This function reads JSON files into whatever type is_type is.
"""
function read_JSON(name_file::String; is_type = Dict{Symbol, Float64})
    nt = nothing
    open(name_file, "r") do f
        nt = JSON2.read(f, is_type)
    end
    nt
end

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

#eventually we will need to go through this and redit load and save
function load_model(file_root::String, p_dict::Dict{Symbol, T}, u_dict::Dict{Symbol, T}; 
            reset_model::Bool = false, gpu::Bool = true, version = :gACh,
            abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
            animate_solution = true, animate_dt = 60.0, 
            model_file_type = :jld2,
            notify = false #set this to true in order to get phone notifications
        ) where T <: Real
    
    if reset_model && isdir(file_root)
        rm(file_root, force = true, recursive = true)
    end

    if !isdir(file_root)
        println("[$(now())]: Making experiment file")
        mkdir(file_root)
    end

    #We first check to see if there is a param loaded
    if isfile("$(file_root)\\params.json")
        p_dict = read_JSON("$(file_root)\\params.json", is_type = Dict{Symbol, T}) 
        p = p_dict |> extract_dict
    else
        #In this scenario we need to use a new dictionary
        println("[$(now())]: Writing a new dictionary")
        p = p_dict |> extract_dict

        write_JSON(p_dict, "$(file_root)\\params.json") #write the parameters to save for later
        write_JSON(u_dict, "$(file_root)\\iconds.json") #write the 
    end
    
    if isfile("$(file_root)\\conds.bson") 
        BSON.@load "$(file_root)\\conds.bson" warmup #This is the BSON file way
        if gpu
            u0 = warmup |> cu
        else
            u0 = warmup
        end
        #Construct the model and ODE problem
        net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ], version = version, gpu = gpu)
        NetProb = SDEProblem(net, noise, u0, (0f0 , 0f1), p)
        
        if model_file_type == :jld2 && !isfile("$(file_root)\\conds.jld2") #In this case we want to make a backup of the file
            JLD2.@save "$(file_root)\\conds.jld2" warmup #This is only here to try to save the older files
        end
    elseif isfile("$(file_root)\\conds.jld2")
        JLD2.@load "$(file_root)\\conds.jld2" warmup #Load the solution from the JSON file
        if gpu
            u0 = warmup |> cu
        else
            u0 = warmup
        end
        #Construct the model and ODE problem
        net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ], version = version, gpu = gpu)
        NetProb = SDEProblem(net, noise, u0, (0f0 , 0f1), p)

        if model_file_type == :bson && !isfile("$(file_root)\\conds.bson") #In this case we want to make a backup of the file
            BSON.@save "$(file_root)\\conds.bson" warmup #as a BSON file
        end 
    else
        println("[$(now())]: Model loaded from new")
        #Load the initial conditions
        if gpu
            u0 = extract_dict(u_dict, p_dict[:nx], p_dict[:ny]) |> cu
        else
            u0 = extract_dict(u_dict, p_dict[:nx], p_dict[:ny])
        end

        #Load the model and ODE interface
        net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ], version = version, gpu = gpu)
        NetProb = SDEProblem(net, noise, u0, (0f0 , p_dict[:t_warm]), p)

        #Warmup the problem
        print("[$(now())]: Warming up the solution... ")
        @time sol = solve(NetProb, SOSRI(), 
                abstol = 2e-2, reltol = 2e-2, maxiters = 1e7,
                progress = true, progress_steps = 1, 
                save_everystep = false
            )
        
        #Save the warmed up solution
        warmup = sol[end]|>Array #Keep this one as the original 
        if model_file_type == :bson
            BSON.@save "$(file_root)\\conds.bson" warmup #This is the BSON file way
        elseif model_file_type == :jld2
            JLD2.@save "$(file_root)\\conds.jld2" warmup #This is only here to try to save the older files
        end

        warmup = sol[end] #Convert this one to CPU
        if gpu
            sol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
        end
    end
    if notify
        BotNotify("{Waves} Model warmup solution loaded")
    end
    
    #Now we want to test whether or not we have run the simulation
    if isfile("$(file_root)\\sol.bson") #Always try to load BSON first
        BSON.@load "$(file_root)\\sol.bson" sol #This is the BSON file way
        if model_file_type == :jld2 && !isfile("$(file_root)\\sol.jld2") #This means we have to backup the file
            JLD2.@save "$(file_root)\\sol.jld2" sol #This is the BSON file way   
        end
    elseif isfile("$(file_root)\\sol.jld2")
        #Load the solution from the JSON file
        JLD2.@load "$(file_root)\\sol.jld2" sol #This is only here to try to save the older files
        if model_file_type == :bson && !isfile("$(file_root)\\sol.bson")
            BSON.@save "$(file_root)\\sol.bson" sol #This is the BSON file way            
        end
    else
        print("[$(now())]: Running the model... ")
        NetProb = SDEProblem(net, noise, warmup, (0f0 , p_dict[:t_run]), p)
        #Run the solution to fruition
        @time sol = solve(NetProb, SOSRI(), 
            abstol = abstol, reltol = reltol, maxiters = maxiters,
            save_idxs = [1:(Int64(p_dict[:nx]*p_dict[:ny]))...], 
            progress = true, progress_steps = 1
        )

        print("[$(now())]: Saving the simulation...")
        sol = convert_to_cpu(sol) #Before saving we need to bump the file to CPU
        if model_file_type == :bson
            BSON.@save "$(file_root)\\sol.bson" sol #This is the BSON file way
        elseif model_file_type == :jld2
            JLD2.@save "$(file_root)\\sol.jld2" sol #This is only here to try to save the older files
        end
        println("Completed")

    end
    if notify
        BotNotify("{Waves} Model solution loaded")
    end

    if !isfile("$(file_root)\\animation.gif") && animate_solution
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

function save_solution(sol, save_path::String; mode = :bson)
    if mode == :bson
        bson("$(save_path)\\sol_data.bson", 
            Dict(
                :sol_t => sol.t, 
                )
        )
        try #This will fail if the file is too big
            
            bson("$(save_path)\\sol_array.bson", 
                Dict(:sol_u => sol.u)
            )
        catch
            println("cannot save something this large. Breaking it up")
            n_half = round(Int64, size(sol.u, 1)/2)
            bson("$(save_path)\\sol_array_pt1.bson", Dict(:sol_u => sol.u[1:n_half-1]))
            bson("$(save_path)\\sol_array_pt2.bson", Dict(:sol_u => sol.u[n_half:end]))
        end
    else
        println("TODO implement JLD2")
    end
end

function load_solution(load_path)
    sol_t = BSON.load("$(load_path)\\sol_data.bson")[:sol_t]
    sol_u = try 
        BSON.load("$(load_path)\\sol_array.bson")[:sol_u]
    catch
        println("Loading method 2")
        sol_array_pt1 = BSON.load("$(load_path)\\sol_array_pt1.bson")[:sol_u]
        sol_array_pt2 = BSON.load("$(load_path)\\sol_array_pt2.bson")[:sol_u]
        vcat(sol_array_pt1, sol_array_pt2)
    end
    p_dict = read_JSON("$(load_path)\\params.json", is_type = Dict{Symbol, Float32})
    p0 = p_dict |> extract_dict
    net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ])
    sol_prob = SDEProblem(net, noise, sol_u[1], (0f0 , p_dict[:t_warm]), p0)
    SciMLBase.build_solution(sol_prob, SOSRI(), sol_t, sol_u) #we can use this to build a solution without GPU
end