"""
This function is used to calculate the current. It represents the generalized form: 

\$I_n = -g_n * R * (V_t - E_n)\$
\$n = (Leak, Ca, K, TREK, ACh)\$
- INPUTS: 
    - Maximal conductance (\$g_n\$) 
    - Gating fraction (\$R_\\infty\$)
    - Membrane Potential (\$V_t\$) 
    - Reversal Potential (\$E_n\$)
- OUTPUTS:
    - Ionic current of type n (\$I_n\$)
"""
I_n(g_n::T, R::T, v::T, E_n::T) where T <: Real = -g_n*R*(v-E_n)

"""
This is the Boltzmann equation for channel gating.  
Origingally described in Hodgkin & Huxley 1955, then simplified into a hyperbolic form by Morris & Lecar 1980

The equation for this is as follows: 

\$ R\\_{\\infty} = \\frac{1}{2} * (1 + tanh(\\frac{V\\_t - V\\_S}{V\\_H}))\$
- INPUTS: 
    - The membrane voltage (\$V_t\$)
    - The slope of the function (\$V_S\$)
    - The half activation point of the function (\$V_H\$)
- OUTPUT: 
    - The gating variable (\$R_{\\infty}\$)

- For Calcium channel activation: \$V_1 = V_S,  V_2 = V_H,  M_{\\infty} = R_{\\infty}\$
- For Potassium channel activation: \$V_3 = V_S,  V_4 = V_H,  N_{\\infty} = R_{\\infty}\$
"""
R_INF(v::T, VS::T, VH::T) where T <: Real  = (1 + tanh((v - VS)/VH))/2;
CuArrays.@cufunc R_INF(v, VS, VH) = (1 + tanh((v - VS)/VH))/2;
R_INF(v::CuArray, VS, VH) = R_INF.(v, VS, VH)

"""
This equation related voltage to the rate constant of opening Potassium channels. Described more in detail in Morris Et. al. 1981. 

\$\\Lambda(V_t, V_3, V_4) =  cosh\\left(\\frac{V_t-V_3}{2V_4} \\right)\$

- INPUTS
    - The membrane voltage (\$V_t\$)
    - The slope value (\$V_3\$)
    - The half saturation value (\$V_4\$)
- OUTPUTS
    - The potassium rate constant (\$\\Lambda\$)
"""
Λ(v::T, V3::T, V4::T) where T <: Real = cosh((v-V3)/(2*V4));
CuArrays.@cufunc Λ(v, V3, V4) = cosh((v-V3)/(2*V4));
Λ(v::CuArray, V3, V4) = LAM.(v, V3, V4) 


Φ(v::T, κ::T, V_0::T) where T <: Real = 1/(1 + exp(-κ * (v - V_0)))
CuArrays.@cufunc Φ(v, κ, V_0) = 1/(1 + exp(-κ * (v - V_0)))
Φ(v::CuArray, κ, V_0) = Φ.(v, κ, V_0)

ħ(e::T, K_d::T) where T <: Real = (e^2)/(e^2 + K_d)
CuArrays.@cufunc ħ(e, K_d) = (e^2)/(e^2 + K_d)
ħ(e::CuArray, K_d) = ħ.(e, K_d)

tar_pars = [:g_leak, :E_leak, :g_Ca, :V1, :V2, :E_Ca, :g_K, :E_K, :g_TREK, :g_ACh, :k_d, :E_ACh, :I_app, :C_m, :V3, :V4, :τn, :C_0, :λ, :δ, :τc, :α, :τa, :β, :τb, :ρ, :τACh, :k, :V0, :σ, :D, :τw]
tar_conds = [:v, :n, :c, :a, :b, :e, :W]
lansdell_pars = [:I_app, :E_Ca, :E_K, :E_Leak, :E_ACh, :V1, :V2, :V3, :V4, :g_Ca, :g_K, :g_Leak, :λ, :g_ACh, :δ, :C_m, :τr, :τs, :τACh, :γ, :α, :β, :k, :V0, :D, :μ]
lansdell_conds = [:v, :r, :s, :a, :W]

#This function is used to get each time step within the simulation
get_timesteps(sol) = [sol.t[i] - sol.t[i-1] for i = 2:length(sol.t)]

include("pde_models.jl")
include("1D_models.jl") 
include("2D_models.jl")    



#Noise models
noise(du::Array{T,1}, u::Array{T,1}, p::Array{T,1}, t::T) where T <: Real = du[end] = p[30]
noise(du::Array{T,3}, u::Array{T,3}, p::Array{T,1}, t::T) where T <: Real = du[:,:,end] .= p[30]
lansdell_noise(du, u, p, t) = du[:,:,end] .= p[end]
