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
u_find(cond::Symbol; list_u::Array{Symbol, 1} = sym_cs) = findall(c -> c == cond, list_u)[1]
#get the index of the parameter in the list
p_find(par::Symbol; list_p::Array{Symbol, 1}  = sym_ps) = findall(p -> p == par, list_p)[1]

"""
When using the Modeling Toolkit, the dictionary needs to be converted into an array of operations
"""
function extract_dict(dict_item::Dict{Symbol, Float64})
    par_set = nothing
    idx = 1
    for key in keys(dict_item)
        if par_set == nothing
            par_set = [Variable(key) => dict_item[key]]
        else
            push!(par_set, (Variable(key)=> dict_item[key]))
        end
    end
    par_set
end

"""
Multiple dispatch function that is called to extract a distribution dictionary
"""
function extract_dict(dict_item::Dict{Symbol, Tuple}, pars::Array{Symbol})
    d0 = Array{Distribution}([])
    for par in pars
        push!(d0, convert_dist(dict_item[par]))
    end
    d0
end

"""
This function walks through the directory and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function parse_abf(super_folder::String; extension::String = ".abf", verbose = false)
    file_list = []
    for (root, dirs, files) in walkdir(super_folder)
        for file in files
            if file[end-3:end] == extension
                path = joinpath(root, file)
                if verbose 
                    println(path) # path to files
                end
                push!(file_list, path)
            end
        end
    end
    file_list
end

function extract_abf(abf_path; verbose = false, v_offset = -25.0)
    if length(abf_path |> splitpath) > 1
        full_path = abf_path
    else
        full_path = joinpath(pwd(), abf_path)   
    end
    #extract the abf file by using pyABF
    exp_data = pyABF.ABF(full_path)
    #if the data is segmented into sweeps (which Jordans data is) concatenate all sweeps
    if length(exp_data.sweepList) > 1
        data = Float64[]
        time = Float64[]
        previous_time = 0.0
        for sweepNumber in exp_data.sweepList
            exp_data.setSweep(sweepNumber = sweepNumber, channel = 0);
            push!(data, exp_data.sweepY...);
            push!(time, (exp_data.sweepX.+previous_time)...);
            previous_time = time[end]
        end
        dt = time[2]*1000
    else
        exp_data.setSweep(sweepNumber = 0, channel = 0);
        data = Float64.(exp_data.sweepY);
        time = Float64.(exp_data.sweepX);
        dt = time[2]*1000
        
    end
    if verbose
        println("Data extracted from $full_path")
        println("Data from time stamp $(t[1]) s to $(t[end]+dt) s with dt = $dt ms")
        println("Data was acquired at $(1/dt/1000) Hz")
        println("$(length(t)) data points")
    end
    data, time, dt
end
