#Used for loading Phys data
using PyCall

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

function load_model(file_root::String, p_dict::Dict{Symbol, T}, u_dict::Dict{Symbol, T}; 
            reset_model::Bool = false, gpu::Bool = true, version = :gACh,
            abstol = 2e-2, reltol = 0.2, maxiters = 1e7,
            animate_solution = true, animate_dt = 60.0, 
            notify = false #set this to true in order to get phone notifications
        ) where T <: Real
    
    if reset_model && isdir(file_root)
        rm(file_root, force = true, recursive = true)
    end

    if !isdir(file_root)
        println("[$(now())]: Making experiment file")
        mkdir(file_root)
    end

    #We can step through each section checking to see what level we are at
    params_file = "$(file_root)\\params.json"
    iconds_file = "$(file_root)\\iconds.json"
    conds_file = "$(file_root)\\conds.bson"
    sol_file = "$(file_root)\\sol.bson"
    animation_file =  "$(file_root)\\animation.gif"

    #We first check to see if there is a param loaded
    if isfile(params_file)
        p_dict = read_JSON(params_file, is_type = Dict{Symbol, T}) 
        p = p_dict |> extract_dict
    else
        #In this scenario we need to use a new dictionary
        println("[$(now())]: Writing a new dictionary")
        p = p_dict |> extract_dict

        write_JSON(p_dict, params_file) #write the parameters to save for later
        write_JSON(u_dict, iconds_file) #write the 
    end
    
    if isfile(conds_file)
        #Load the solution from the JSON file
        BSON.@load conds_file warmup #This is the BSON file way
        u0 = warmup |> cu
        #Construct the model and ODE problem
        net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ], version = version, gpu = gpu)
        NetProb = SDEProblem(net, noise, u0, (0f0 , 0f1), p)
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
        BSON.@save conds_file warmup #as a JSON
        warmup = sol[end] #Convert this one to CPU
        if gpu
            sol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
        end
    end
    if notify
        BotNotify("{Waves} Model warmup solution loaded")
    end
    
    #Now we want to test whether or not we have run the simulation
    if isfile(sol_file)
        #We need to have previously activated the problem
        BSON.@load sol_file sol #This is the JSON method
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
        BSON.@save sol_file sol #Save as a BSON
        println("Completed")

    end
    if notify
        BotNotify("{Waves} Model solution loaded")
    end

    if !isfile(animation_file) && animate_solution
        println("[$(now())]: Animating simulation...")
        anim = @animate for t = 1.0:animate_dt:sol.t[end]
            frame_i = reshape(sol(t) |> Array, (p_dict[:nx]|>Int64, p_dict[:ny]|>Int64))
            heatmap(frame_i, ratio = :equal, grid = false,
                    xaxis = "", yaxis = "", xlims = (0, p_dict[:nx]), ylims = (0, p_dict[:ny]),
                    c = :curl, clims = (-70.0, 0.0),
            )
        end
        gif(anim, animation_file, fps = 1000.0/animate_dt)
    end

    #We want to return the solution
    return sol
end 