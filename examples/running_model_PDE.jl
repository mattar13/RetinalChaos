using Plots
using DifferentialEquations, ModelingToolkit
using MethodOfLines, DomainSets, DiffEqOperators

@parameters x y t
@variables v̂(..) n̂(..) m̂(..) ĥ(..) ĉ(..) â(..) b̂(..) ê(..) î(..) Ŵ(..)
@parameters g_leak E_leak g_Ca V1 V2 E_Ca g_K E_K g_TREK g_ACh k_ACh E_ACh g_GABA k_GABA E_Cl I_app C_m
@parameters V3 V4 τn
@parameters C_0 λ δ τc
@parameters α τa β τb ρe ρi τACh τGABA VSe VSi V0e V0i
@parameters De Di
@parameters τw σ
@parameters g_Na E_Na
@parameters V7 V8 V9
@parameters V10 V11 V12
@parameters V13 V14 V15
@parameters V16 V17 V18

Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

println("Loading auxillary functions")
M∞(v̂) = (1 + tanh((v̂ - V1) / V2)) / 2
N∞(v̂) = (1 + tanh((v̂ - V3) / V4)) / 2
Λ(v̂) = cosh((v̂ - V3) / (2 * V4));
Φe(v̂) = 1 / (1 + exp(-VSe * (v̂ - V0e)))
Φi(v̂) = 1 / (1 + exp(-VSi * (v̂ - V0i)))
ħe(ê) = (ê^2) / (ê^2 + k_ACh)
ħi(î) = (î^2) / (î^2 + k_GABA)

α_M(v̂) = -(v̂ - V8) / (V7 * (exp(-(v̂ - V8) / V9) - 1))
β_M(v̂) = V10 * (exp(-(v̂ - V11) / V12))
α_H(v̂) = V13 * (exp(-(v̂ - V14) / V15))
β_H(v̂) = 1 / (V16 * (exp(-(v̂ - V17) / V18) + 1))

∇²(u) = Dxx(u) + Dyy(u) #This is the diffusion aspect of the equation

#Each one of these versions needs a PDE version
ÎLeak(v̂) = -g_leak * (v̂ - E_leak)
ÎCa(v̂) = -g_Ca * M∞(v̂) * (v̂ - E_Ca)
ÎK(v̂) = -g_K * n̂(x, y, t) * (v̂ - E_K)
ÎTREK(v̂) = -g_TREK * b̂(x, y, t) * (v̂ - E_K)
ÎACh(v̂) = -g_ACh * ħe(ê(x, y, t)) * (v̂ - E_ACh)
ÎGABA(v̂) = -g_GABA * ħi(î(x, y, t)) * (v̂ - E_Cl)
ÎNa(v̂) = -g_Na * m̂(x, y, t)^3 * ĥ(x, y, t) * (v̂ - E_Na)

PDEeqs = [
     Dt(v̂(x, y, t)) ~ (ÎLeak(v̂(x, y, t)) + ÎCa(v̂(x, y, t)) + ÎK(v̂(x, y, t)) + ÎTREK(v̂(x, y, t)) + ÎACh(v̂(x, y, t)) + ÎGABA(v̂(x, y, t)) + ÎNa(v̂(x, y, t)) + I_app + Ŵ(x, y, t)) / C_m,
     Dt(n̂(x, y, t)) ~ (Λ(v̂(x, y, t)) * ((N∞(v̂(x, y, t)) - n̂(x, y, t)))) / τn,
     Dt(m̂(x, y, t)) ~ α_M(v̂(x, y, t)) * (1 - m̂(x, y, t)) - β_M(v̂(x, y, t)) * m̂(x, y, t),
     Dt(ĥ(x, y, t)) ~ α_H(v̂(x, y, t)) * (1 - ĥ(x, y, t)) - β_H(v̂(x, y, t)) * ĥ(x, y, t),
     Dt(ĉ(x, y, t)) ~ (C_0 + δ * ÎCa(v̂(x, y, t)) - λ * ĉ(x, y, t)) / τc,
     Dt(â(x, y, t)) ~ (α * ĉ(x, y, t)^4 * (1 - â(x, y, t)) - â(x, y, t)) / τa,
     Dt(b̂(x, y, t)) ~ (β * â(x, y, t)^4 * (1 - b̂(x, y, t)) - b̂(x, y, t)) / τb,
     Dt(ê(x, y, t)) ~ (De * ∇²(ê(x, y, t)) + ρe * Φe(v̂(x, y, t)) - ê(x, y, t)) / τACh,
     Dt(î(x, y, t)) ~ (Di * ∇²(î(x, y, t)) + ρi * Φi(v̂(x, y, t)) - î(x, y, t)) / τGABA,
     Dt(Ŵ(x, y, t)) ~ -Ŵ(x, y, t) / τw
]

#set the time parameters 
dx = 0.5
dy = 0.5
dt = 1.0
xmin = ymin = tmin = 0.0
xmax = ymax = 3.0
tmax = 120e3
tsteps = tmin:dt:tmax

ics = Dict(
     :v => -65.0,
     :n => 0.0,
     :c => 0.085,
     :a => 0.1,
     :b => 0.03,
     :e => 0.0,
     :i => 0.0,
     :W => 0.0
)

#Load the boundary conditions
bcs = [
     #Time boundary conditions
     v̂(x, y, tmin) ~ -65.0,
     n̂(x, y, tmin) ~ 0.0,
     m̂(x, y, tmin) ~ 0.0,
     ĥ(x, y, tmin) ~ 0.0,
     ĉ(x, y, tmin) ~ 0.085,
     â(x, y, tmin) ~ 0.1,
     b̂(x, y, tmin) ~ 0.03,
     ê(x, y, tmin) ~ 0.0,
     î(x, y, tmin) ~ 0.0,
     Ŵ(x, y, tmin) ~ 0.0,
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

#Set the domains
domains = [
     t ∈ Interval(tmin, tmax)
     x ∈ Interval(xmin, xmax)
     y ∈ Interval(ymin, ymax)
]

ivs = [x, y, t]
dvs = [v̂(x, y, t), n̂(x, y, t), m̂(x, y, t), ĥ(x, y, t), ĉ(x, y, t), â(x, y, t), b̂(x, y, t), ê(x, y, t), î(x, y, t), Ŵ(x, y, t)]

parameters = [
     g_leak => 2.0,
     E_leak => -70.0,
     g_Ca => 8.5,
     V1 => -20.0,
     V2 => 20.0,
     E_Ca => 50.0,
     g_K => 4.0,
     E_K => -90.0,
     g_TREK => 3.0,
     I_app => 0.0,
     C_m => 13.6,
     V3 => -25.0,
     V4 => 7.0,
     τn => 5.0,
     C_0 => 0.088,
     λ => 2.702,
     δ => 0.010503,
     τc => 2000.0,
     α => 625.0,
     τa => 8300.0,
     β => 34.0,
     τb => 8300.0,
     ρe => 6.0,
     ρi => 5.0,
     τACh => 540.0,
     τGABA => 1000.0,
     VSe => 0.2,
     VSi => 0.2,
     V0e => -40.0,
     V0i => -40.0,
     g_ACh => 0.215,
     k_ACh => 0.1,
     E_ACh => 0.0,
     g_GABA => 0.9,
     k_GABA => 0.1,
     E_Cl => -55.0,
     g_Na => 1.0,
     E_Na => 55.0, #55.0
     De => 0.005,
     Di => 0.005,
     τw => 800.0,
     σ => 0.1,
     V7 => 10.0,
     V8 => -40.0,
     V9 => 10.0,
     V10 => 4.0,
     V11 => -65.0,
     V12 => 18.0,
     V13 => 0.07,
     V14 => -65.0,
     V15 => 20.0,
     V16 => 1.0,
     V17 => -35.0,
     V18 => 10.0,
]

@named SAC_PDE = PDESystem(PDEeqs, bcs, domains, ivs, dvs, parameters) #Create the undiscretized PDE system
discretization = MOLFiniteDifference([x => dx, y => dy], t)
grid = get_discrete(SAC_PDE, discretization) #Make a representation of the discrete map
discrete_x = grid[x]
discrete_y = grid[y]
nx = length(discrete_x)
ny = length(discrete_y)
@time probPDE = discretize(SAC_PDE, discretization) # This gives an ODEProblem since it's time-dependent

sig = 0.1
function g(du, u, p, t)
     du[end-(nx*ny):end] .= sig
     return du
end
probSPDE = SDEProblem(probPDE.f, g, probPDE.u0, probPDE.tspan, probPDE.p)

#%% Run the model
@time solSDE = solve(probSPDE, SOSRI(), progress=true, progress_steps=1)

#%% Plot all the answers
discrete_t = solSDE.t

anim = @animate for i in 1:25:length(discrete_t)
     println(i)
     v_map = map(d -> solSDE[d][i], grid[v(t, x, y)])
     h1 = heatmap(discrete_x, discrete_y, v_map, clim=(-90.0, 10.0), title="$(discrete_t[i])", aspect_ratio=:equal)
end
gif(anim, "sync_voltage.gif", fps=60.0)

x_traces = reshape(solSDE[grid[v̂(x, y, t)]], nx*ny)
Plots.plot(discrete_t, x_traces)