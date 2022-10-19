function data_bytesize(arr::Vector{T}; bson_header_bytes = 133) where T <: Real
    data_bytes = sizeof(arr)
    size_bytes = 11 #There is a section in the header that 
    file_size = bson_header_bytes + data_bytes + size_bytes
    return file_size
end

function data_bytesize(arr::Matrix{T}; bson_header_bytes = 133) where T <: Real
    data_bytes = sizeof(arr)
    size_bytes = 11 * ndims(arr) #There is a section in the header that 
    file_size = bson_header_bytes + data_bytes + size_bytes
    return file_size
end

function data_bytesize(arr::Vector{Vector{T}}; bson_header_bytes = 133) where T <: Real
    data_bytes = sizeof(arr[1])
    size_bytes = 111 * length(arr)
    return bson_header_bytes + data_bytes + 111*length(arr)
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
> Out[2]: 125x125xN AbstractArray:
"""
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}) = map(x -> Float64(dict_item[x]), pars)
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}, dims::Tuple) = cat(map(x -> fill(Float64(dict_item[x]), dims), pars)..., dims = length(dims)+1)

#get the index of the conditions in the list
indexof(num::Symbolics.Num; syms = ODEModel.ps) = findfirst(isequal(num), syms)
indexof(sym::Symbol; syms = ODEModel.ps) = findfirst(isequal(ODEModel.var_to_name[sym]), syms)
"""
This function can do a few things

    1) It can take a dictionary and extract it in the order of the dict_names
    2) It can read through a dictionary and search for certain keys
"""
function extract_dict(dict_item::Dict{Symbol,T}, key_names; dims = (1)) where {T<:Real}
    #This is a parameter and we should extract it
    items = T[]
    for key in key_names
        push!(items, dict_item[key])
    end
    if dims == 1
        return items
    else
        val_map = zeros(dims..., length(items))
        for i in 1:length(items)
            val_map[:, :, i] .= items[i]
        end
        return val_map
    end
end

function extract_dict(dict_item::Dict{Symbol, T}; key_library = [t_pars, t_conds]) where T <: Real
    #To extract the dictionary first we want to take a look at the first param to see if it is a parameter or a condition
    for key_names in key_library
        #println("Here")
        if dict_item.keys[1] ∈ key_names
            items = extract_dict(dict_item, key_names)
            return items
        end
    end
    println("if you made it here there was an error")
    throw("NoKeyLibraries")
end

function extract_dict(dict_item::Dict{Symbol, T}, dims...; dict_names = [t_pars, t_conds]) where T <: Real
    #Extract the dict_item
    vals = extract_dict(dict_item; dict_names = dict_names)
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
function save_solution(sol_t::Vector{T}, sol_u::Array{Vector{T}, 1}, save_path::String; 
    name = "sol", partitions = 1, mode = :bson
    ) where T <: Real
    
    if mode == :bson
        if partitions == -1 #This means data is getting appended
            #check to see if a file_contents file exists
            file_contents = read_JSON(Dict{Symbol, Vector{String}}, "$(save_path)\\file_contents.json")
            push!(file_contents[:t], "$(name)_t.bson")
            push!(file_contents[:u], "$(name)_u.bson")
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            bson("$(save_path)\\$(name)_t.bson", Dict(:t => sol_t))
            bson("$(save_path)\\$(name)_u.bson", Dict(:u => sol_u))
        elseif partitions == 1
            file_contents = Dict(:t => ["$(name)_t.bson"], :u => ["$(name)_u.bson"])
            #write the details to a 
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            bson("$(save_path)\\$(name)_t.bson", Dict(:t => sol_t))
            bson("$(save_path)\\$(name)_u.bson", Dict(:u => sol_u))
        else
            file_contents = Dict(:t => ["$(name)1_t.bson"], :u => ["$(name)1_u.bson"])
            for i in 2:partitions
                push!(file_contents[:t], "$(name)$(i)_t.bson")
                push!(file_contents[:u], "$(name)$(i)_u.bson")
            end
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            partition_idxs = round.(Int64, LinRange(1, length(sol_t), partitions+1))
            for idx in 1:length(partition_idxs)-1
                if idx == 1
                    bson("$(save_path)\\$(name)$(idx)_t.bson", Dict(:t => sol_t[1:partition_idxs[idx+1]]))
                    bson("$(save_path)\\$(name)$(idx)_u.bson", Dict(:u => sol_u[1:partition_idxs[idx+1]]))
                else
                    bson("$(save_path)\\$(name)$(idx)_t.bson", Dict(:t => sol_t[partition_idxs[idx]+1:partition_idxs[idx+1]]))
                    bson("$(save_path)\\$(name)$(idx)_u.bson", Dict(:u => sol_u[partition_idxs[idx]+1:partition_idxs[idx+1]]))
                end
            end
        end
    else
        println("TODO implement JLD2")
    end
end

save_solution(sol, save_path::String; kwargs...) = save_solution(sol.t, sol.u; karwgs...)
save_solution(res::DiffEqCallbacks.SavedValues{Float64, Vector{Float64}}, save_path::String; kwargs...) = save_solution(res.t, res.saveval, save_path; kwargs...)

function load_solution(load_path)
    file_contents = read_JSON(Dict{Symbol, Vector{String}}, "$(load_path)\\file_contents.json")
    n_partitions = length(file_contents[:u])
    if n_partitions == 1 #there are no partitions in the file
        sol_t = BSON.load("$(load_path)\\$(file_contents[:t][1])")[:t]
        sol_u = BSON.load("$(load_path)\\$(file_contents[:u][1])")[:u]
    else
        sol_t = Float64[]
        sol_u = Vector{Vector{Float64}}()
        for i in 1:n_partitions
            sol_ti = BSON.load("$(load_path)\\$(file_contents[:t][i])")[:t]
            sol_ui = BSON.load("$(load_path)\\$(file_contents[:u][i])")[:u]
            push!(sol_t, sol_ti...)
            push!(sol_u, sol_ui...)
        end
    end

    p_dict = read_JSON(Dict{Symbol, Float32}, "$(load_path)\\params.json")
    p0 = p_dict |> extract_dict
    net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ])
    sol_prob = SDEProblem(net, noise, sol_u[1], (0f0 , p_dict[:t_warm]), p0)
    SciMLBase.build_solution(sol_prob, SOSRI(), sol_t, sol_u) #we can use this to build a solution without GPU
end

"""
This function runs the model using the indicated parameters

"""
function run_model(file_root::String, p_dict::Dict{Symbol, T}, u_dict::Dict{Symbol, T}; 
        gpu::Bool = true, version::Symbol = :gACh,
        abstol::Float64 = 2e-2, reltol::Float64 = 2e-2, maxiters::Float64 = 1e7,
        model_file_type = :bson, 
        iterations = 1 #this option sets the model up into sections so that we can break up the saving of the solutions
    ) where T <: Real
    
    if !isdir(file_root)
        mkdir(file_root)
    end

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
    #Add this in because it seems that you need to constantly clean the function
    sol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()

    #Save the warmed up solution
    if model_file_type == :bson #In this case we want to make a backup of the file
        BSON.@save "$(file_root)\\conds.bson" warmup #as a BSON file
    elseif model_file_type == :jld2  #In this case we want to make a backup of the file
        JLD2.@save "$(file_root)\\conds.jld2" warmup #This is only here to try to save the older files
    end

    results = RetinalChaos.SavedValues(Float64, Vector{Float64})
    cb = RetinalChaos.SavingCallback(
        (u, t, integrator) -> reshape(u[:,:,1], size(u,1) * size(u,2)), 
        results
    )

    print("[$(now())]: Running the model... ")
    NetProb = SDEProblem(net, noise, warmup, (0f0 , p_dict[:t_run]), p0)
    #Run the solution saving values to results
    @time sol = solve(NetProb, SOSRI(), 
        callback = cb, #This saves the solution without actually saving anything to the GPU
        abstol = abstol, reltol = reltol, maxiters = maxiters,
        save_everystep = false, 
        progress = true, progress_steps = 1
    )
    #make sure to zero out the solution to save GPU space 
    sol = nothing; GC.gc(true); RetinalChaos.CUDA.reclaim()
    
    #We want to return the solution
    return results
end 