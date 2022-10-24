"""
This function returns a ODE function ready to go
"""
function loadODE(u0::Dict{Num,T}, parameters::Dict{Num,T}; tmin = 0.0, tmax=300e3) where {T<:Real}
     return ODEProblem(ODEModel, u0, (tmin, tmax), parameters)
end

function loadSDE(u0::Dict{Num,T}, parameters::Dict{Num,T}; tmin = 0.0, tmax=300e3) where {T<:Real}
     return SDEProblem(SDEModel, u0, (tmin, tmax), parameters)
end

function loadPDE(u0::Dict{Num,T}, parameters::Dict{Num,T};
     tmin=0.0, dt=1.0, tmax=300e3,
     xmin=0.0, dx=0.5, xmax=3.0,
     ymin=0.0, dy=0.5, ymax=3.0
) where {T<:Real}
     print("Setting up PDE model... ")
     bcs = [
          #Time boundary conditions
          v̂(x, y, tmin) ~ u0[v],
          n̂(x, y, tmin) ~ u0[n],
          m̂(x, y, tmin) ~ u0[m],
          ĥ(x, y, tmin) ~ u0[h],
          ĉ(x, y, tmin) ~ u0[c],
          â(x, y, tmin) ~ u0[a],
          b̂(x, y, tmin) ~ u0[b],
          ê(x, y, tmin) ~ u0[e],
          î(x, y, tmin) ~ u0[i],
          Ŵ(x, y, tmin) ~ u0[W],
          #Spatial Boundary Conditions for the variable e
          Dx(ê(xmin, y, t)) ~ 0.0,
          Dx(ê(xmax, y, t)) ~ 0.0,
          Dy(ê(x, ymin, t)) ~ 0.0,
          Dy(ê(x, ymax, t)) ~ 0.0,
          Dx(î(xmin, y, t)) ~ 0.0,
          Dx(î(xmax, y, t)) ~ 0.0,
          Dy(î(x, ymin, t)) ~ 0.0,
          Dy(î(x, ymax, t)) ~ 0.0,
     ]
     domains = [
          x ∈ Interval(xmin, xmax)
          y ∈ Interval(ymin, ymax)
          t ∈ Interval(tmin, tmax)
     ]
     dimensions = [x, y, t]
     states = [v̂(x, y, t), n̂(x, y, t), m̂(x, y, t), ĥ(x, y, t), ĉ(x, y, t), â(x, y, t), b̂(x, y, t), ê(x, y, t), î(x, y, t), Ŵ(x, y, t)]
     ps = Vector{Pair{Num,Float64}}()
     for entry in parameters
          push!(ps, entry[1] => entry[2])
     end
     @named PDEModel = PDESystem(PDEeqs, bcs, domains, dimensions, states, ps) #Create the undiscretized PDE system
     println("Complete")
     print("Discretizing model... ")
     discretization = MOLFiniteDifference([x => dx, y => dy], t)
     # This gives an ODEProblem since it's time-dependent
     grid = get_discrete(PDEModel, discretization) #Make a representation of the discrete map
     discreteX = grid[x]
     discreteY = grid[y]
     discreteT = tmin:dt:tmax
     #Discretize the ODE
     @time probPDE = discretize(PDEModel, discretization)
     println("Complete")
     return grid, probPDE
end

function loadSPDE(u0::Dict{Num,T}, parameters::Dict{Num,T}; kwargs...) where {T<:Real}
     grid, probPDE = loadPDE(u0, parameters; kwargs...)
     print("Adding a noise term to the equation... ")
     sig = 0.1
     function g(du, u, p, t)
          du[end-(nx*ny):end] .= sig
          return du
     end
     probSPDE = SDEProblem(probPDE.f, g, probPDE.u0, probPDE.tspan, probPDE.p)
     println("Complete")
     return grid, probSPDE
end