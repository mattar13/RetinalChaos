#using Revise
#using RetinalChaos

#%% Test out diffusion at different biases
using Plots
using DifferentialEquations, ModelingToolkit
using MethodOfLines, DomainSets, DiffEqOperators


xmin = 0.0
ymin = 0.0
tmin = 0.0
dx = 1.0
dy = 1.0
dt = 1.0
xmax = 10.0
ymax = 10.0
tmax = 400.0
tstops = tmin:dt:tmax

@parameters t x y
@variables e(..)
@parameters De

function ∇²(u; l = 2.0, r = 2.0, up = 1.0, down = 1.0)
     Dyy(u) + Dxx(u) #This is the diffusion aspect of the equation
end

Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

PDEeqs = [
     Dt(e(x, y, t)) ~ De * ∇²(e(x, y, t)) #- e(x,y,t)/100.0,
     #Dt(e(x, t)) ~ 0.1 * ∇²(e(x, t))
]

bcs = [
     #Time boundary conditions
     e(x, y, tmin) ~ 0.0,
     #Spatial Boundary Conditions for the variable e
     Dx(e(xmin, y, t)) ~ 0.0,
     Dx(e(xmax, y, t)) ~ 0.0,
     #Dx(e(xmin, t)) ~ 0.0,
     #Dx(e(xmax, t)) ~ 0.0,
     Dy(e(x, ymin, t)) ~ 0.0,
     Dy(e(x, ymax, t)) ~ 0.0,
]

domains = [
     x ∈ Interval(xmin, xmax)
     y ∈ Interval(ymin, ymax)
     t ∈ Interval(tmin, tmax)
]

dimensions = [x, y, t]
states = [e(x, y, t)]

ps = [De => 0.01]

@named PDEModel = PDESystem(PDEeqs, bcs, domains, dimensions, states, ps) #Create the undiscretized PDE system
discretization = MOLFiniteDifference([x => dx, y => dy], t)
@time probPDE = discretize(PDEModel, discretization)
grid = get_discrete(PDEModel, discretization) #Make a representation of the discrete map

# Discretization complete
new_u0 = probPDE.u0
size(probPDE.u0)
new_u0[50] = 1.0
new_prob = remake(probPDE; u0=new_u0)
# Run the model
sol = solve(new_prob, saveat=tstops, progress=true, progress_steps=1)

#%%
discrete_t = sol.t
discrete_x = grid[x]
discrete_y = grid[y]
anim = @animate for i in 1:10:length(discrete_t)
     println(i)
     v_map = map(d -> sol[d][i], grid[e(x, y, t)])
     h1 = heatmap(discrete_x, discrete_y, v_map, clim=(0.0, 0.6), title="$(discrete_t[i])", aspect_ratio=:equal)
end
gif(anim, "test_diffusion.gif", fps=10.0)

#%%
e_traces = reshape(sol[grid[e(x, y, t)]], length(discrete_x) * length(discrete_y))
p2 = Plots.plot(discrete_t, e_traces)