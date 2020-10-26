#All of our modules from differential equations can be imported
import RetinalChaos: SDEProblem, ODEProblem, solve, SOSRI, @ode_def
import RetinalChaos: write_JSON, read_JSON, extract_dict, Network
import RetinalChaos: lansdell_pars, lansdell_conds, lansdell_noise
using Plots
import DifferentialEquations: NoiseProblem, NoiseFunction

cd("C:\\Users\\mtarc\\JuliaScripts\\RetinalChaos\\Notebooks")
dims = (64,64)
lansdell_net = Network(dims...; version = :Lansdell);
p_dict = read_JSON("lansdell_params.json"; is_type = Dict{Symbol, Float64});
u_dict = read_JSON("lansdell_conds.json"; is_type = Dict{Symbol, Float64});
p0 = extract_dict(p_dict, lansdell_pars);
u0 = extract_dict(u_dict, lansdell_conds, dims);
tspan = (0.0, 100e3)
prob = SDEProblem(lansdell_net, lansdell_noise, u0, tspan, p0);
sol = solve(prob, SOSRI(), abstol = 0.2, reltol = 0.2, maxiter = 1e7, progress = true, saveat = 100.0)
p = plot(layout = grid(4,1), legend = :false)
plot(sol.t./1000, sol[1,1,1,:], ylabel = :Vt)
plot!(p[2], sol.t./1000, sol[1,1,2,:], ylabel = :Rt)
plot!(p[3], sol.t./1000, sol[1,1,3,:], ylabel = :St)
plot!(p[4], sol.t./1000, sol[1,1,4,:], ylabel = :At, xlabel = "Time (ms)")
#plot!(p[5], sol.t./1000, sol[1,1,5,:])
p


import RetinalChaos.frame_draw
frame_draw(Array(sol); idx = 1)



import RetinalChaos.calculate_threshold

import RetinalChaos.fIn
Ui = rand(1000)
Wi = randn(1000)*0.25
Wt = map(t->sol(t)[1,1,5], sol.t[1:10:10000])

Ui_exp = map(w -> fIn(abs(w), 2.0, -60.0, 50.0), Ui)
Wi_exp = map(w -> fIn(abs(w), 2.0, -60.0, 50.0), Wi)
Wt_exp = map(w -> fIn(abs(w), 2.0, -60.0, 50.0), Wt)

p = plot(layout = grid(4, 3), legend = false)
plot!(p[1,1], Ui, ylims = (minimum(Wt),maximum(Wt)))
plot!(p[1,2], Wi, ylims = (minimum(Wt),maximum(Wt)))
plot!(p[1,3], Wt, ylims = (minimum(Wt),maximum(Wt)))

histogram!(p[2,1], Ui)
histogram!(p[2,2], Wi)
histogram!(p[2,3], Wt)

plot!(p[3,1], Ui_exp, ylims = (minimum(Wt_exp),maximum(Wt_exp)))
plot!(p[3,2], Wi_exp, ylims = (minimum(Wt_exp),maximum(Wt_exp)))
plot!(p[3,3], Wt_exp, ylims = (minimum(Wt_exp),maximum(Wt_exp)))

histogram!(p[4,1], Ui_exp)
histogram!(p[4,2], Wi_exp)
histogram!(p[4,3], Wt_exp)
p

##Notes from lansdell
#The kind of noise used is
using Plots
import RetinalChaos: extract_dict, read_JSON
import RetinalChaos: ODEProblem, SDEProblem, SOSRI, solve
import RetinalChaos: TModel_sys
p_dict = read_JSON("params.json", is_type = Dict{Symbol, Float64});
u_dict = read_JSON("conds.json", is_type = Dict{Symbol, Float64});
tspan = (0.0, 300e3)
u0 = extract_dict(u_dict, TModel_sys)
p0 = extract_dict(p_dict, TModel_sys)
prob = SDEProblem(TModel_sys, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), abstol = 0.2, reltol = 0.2, maxiters = 1e7, progress = true)
plot(sol, vars = [:v, :c, :ACh], layout = grid(3, 1))
