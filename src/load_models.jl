"""
This function returns a ODE function ready to go
"""
function loadODE(u0::Dict{Num, T}, parameters::Dict{Num, T}; tmax=300e3) where T <: Real
     tmin = 0.0
     return ODEProblem(ODEModel, u0, (tmin, tmax), parameters)
end

function loadSDE(u0::Dict{Num,T}, parameters::Dict{Num,T}; tmax=300e3) where {T<:Real}
     tmin = 0.0
     return SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
end

function loadPDE(u0::Dict{Num, T}, parameters::Dict{Num, T};
     tmax = 300e3, xmax = 3.0, ymax = 3.0 
) where T<:Real


end