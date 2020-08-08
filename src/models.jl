"""
This is the Boltzmann equation for channel gating. 
Origingally described in Hodgkin & Huxley Et al 1955, then simplified into a hyperbolic form by Morris & Lecar 1980

The equation for this is as follows: 

\$ R\\_{\\infty} = \\frac{1}{2} * (1 + tanh(\\frac{V\\_t - V\\_S}{V\\_H}))\$

For Calcium channel activation we can use: \$M\\_{\\infty}\$

For Potassium channel activation we can use \$N\\_{\\infty}\$

To use this function
"""
R_INF(V, VS, VH) = (1 + tanh((V - VS)/VH))/2;


Λ(V, V3, V4) = cosh((V-V3)/(2*V4));
Φ(v, κ, V_0) = 1/(1 + exp(-κ * (v - V_0)))
ħ(e, K_d) = (e^2)/(e^2 + K_d)
h(e, δ) = (δ*e^2)/(1+δ*e^2)

CuArrays.@cufunc M_INF(V, V1, V2) = (1 + tanh((V - V1)/V2))/2;
CuArrays.@cufunc N_INF(V, V3, V4) = (1 + tanh((V - V3)/V4))/2;
CuArrays.@cufunc Λ(V, V3, V4) = cosh((V-V3)/(2*V4));
CuArrays.@cufunc Φ(v, κ, V_0) = 1/(1 + exp(-κ * (v - V_0)))
CuArrays.@cufunc ħ(e, K_d) = (e^2)/(e^2 + K_d)
CuArrays.@cufunc h(e, δ) = (δ*e^2)/(1+δ*e^2)

M_INF(V::CuArray, V1, V2) = M_INF.(V, V1, V2)
N_INF(V::CuArray, V3, V4) = N_INF.(V, V3, V4)
Λ(V::CuArray, V3, V4) = LAM.(V, V3, V4) 
Φ(v::CuArray, κ, V_0) = Φ.(v, κ, V_0)
ħ(e::CuArray, K_d) = ħ.(e, K_d)
h(e::CuArray, δ) = h.(e, δ)

#This is the generalized current calculation function
fI(g::Float64, r, v, e::Float64) = -g*r*(v-e)
#This is for calculating filtered shot noise
fIn(W::Float64, λ::Float64, v, e::Float64) = -(-log(W)/λ)*(v-e)  


#p_find(p) = findall(x -> x == p, BurstModel.params)
#u_find(u) = findall(x -> x == u, BurstModel.syms)
tar_pars = [:g_leak, :E_leak, :g_Ca, :V1, :V2, :E_Ca, :g_K, :E_K, :g_TREK, :g_ACh, :k_d, :E_ACh, :I_app, :C_m, :V3, :V4, :τn, :C_0, :λ, :δ, :τc, :α, :τa, :β, :τb, :ρ, :τACh, :k, :V0, :σ, :D, :τw]
tar_conds = [:v, :n, :c, :a, :b, :e, :W]
lansdell_pars = [:I_app, :E_Ca, :E_K, :E_Leak, :E_ACh, :V1, :V2, :V3, :V4, :g_Ca, :g_K, :g_Leak, :λ, :g_ACh, :δ, :C_m, :τr, :τs, :τACh, :γ, :α, :β, :k, :V0, :D, :μ]
lansdell_conds = [:v, :r, :s, :a, :W]

#This function is used to get each time step within the simulation
get_timesteps(sol) = [sol.t[i] - sol.t[i-1] for i = 2:length(sol.t)]

include("pde_models.jl")
include("1D_models.jl") 
include("2D_models.jl")    

sym_ps = map(x -> x |> Symbol, T_ode.ps)
sym_cs = map(x -> x |> Symbol, T_ode.states)

@register M_INF(v, V1, V2)
@register N_INF(v, V3, V4)
@register Λ(v, V3, V4)
@register Φ(v, κ, V0)
@register ħ(e, k_d)
@register h(e, δ)

#Noise models
noise(du, u, p, t) = du[end] = p[30]
noise_2D(du, u, p, t) = du[:,:,end] .= p[30]
lansdell_noise(du, u, p, t) = du[:,:,end] .= p[end]

