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
        if isnothing(par_set)
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