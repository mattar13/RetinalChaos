using Plots
using DifferentialEquations, DiffEqOperators
using BandedMatrices

#%% Using centered difference and banded matrices
deriv = 2
order = 2
Δx = 0.001
N = 100
Δ = CenteredDifference(deriv, order, Δx, N)

u0 = zeros(N)
u0[1:5] .= 10.0
u0[6:10] .= 1.0
c = LinRange(10.0, 1.0, N-10) |> collect
c = [fill(10, 5); fill(1.0, 5); c]
bc = Neumann0BC(Δx)
#bc = DirichletBC(100.0, 0.0)
function step(du, u,p,t) 
    du .= min.(Δ*bc*u, c-u)
    du[1:5] .= (c[1:5] - u[1:5]) 
end
prob = ODEProblem(step, u0, (0.0, 1000.0); progress = true)
sol = solve(prob)
#%%
heatmap(sol |> Array)
#%%
t = collect(1:sol.t[end])
#%%
anim = @animate for t in sol.t
    plot(sol(t), ylims = (0, maximum(sol)))
end
gif(anim)