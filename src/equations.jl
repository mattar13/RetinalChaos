println("Loading auxillary functions")
M∞(v) = (1 + tanh((v - V1) / V2)) / 2
N∞(v) = (1 + tanh((v - V3) / V4)) / 2
Λ(V) = cosh((V - V3) / (2 * V4));
Φe(v) = 1 / (1 + exp(-VSe * (v - V0e)))
Φi(v) = 1 / (1 + exp(-VSi * (v - V0i)))
ħe(e) = (e^2) / (e^2 + k_ACh)
ħi(i) = (i^2) / (i^2 + k_GABA)

α_M(v) = -(v - V8) / (V7 * (exp(-(v - V8) / V9) - 1))
β_M(v) = V10 * (exp(-(v - V11) / V12))
α_H(v) = V13 * (exp(-(v - V14) / V15))
β_H(v) = 1 / (V16 * (exp(-(v - V17) / V18) + 1))

ILeak(v) = -g_leak * (v - E_leak)
ICa(v) = -g_Ca * M∞(v) * (v - E_Ca)
IK(v) = -g_K * n * (v - E_K)
ITREK(v) = -g_TREK * b * (v - E_K)
IACh(v) = -g_ACh * ħe(e) * (v - E_ACh)
IGABA(v) = -g_GABA * ħi(i) * (v - E_Cl)
INa(v) = -g_Na * m^3 * h * (v - E_Na)

#bias_yP
#bias_
function ∇²(u)
     println(size(u))
     Dyy(u) - (2 * u) + Dxx(u) #This is the diffusion aspect of the equation
end

#Each one of these versions needs a PDE version
ÎLeak(v̂) = -g_leak * (v̂ - E_leak)
ÎCa(v̂) = -g_Ca * M∞(v̂) * (v̂ - E_Ca)
ÎK(v̂) = -g_K * n̂(x, y, t) * (v̂ - E_K)
ÎTREK(v̂) = -g_TREK * b̂(x, y, t) * (v̂ - E_K)
ÎACh(v̂) = -g_ACh * ħe(ê(x, y, t)) * (v̂ - E_ACh)
ÎGABA(v̂) = -g_GABA * ħi(î(x, y, t)) * (v̂ - E_Cl)
ÎNa(v̂) = -g_Na * m̂(x, y, t)^3 * ĥ(x, y, t) * (v̂ - E_Na)

# These equations are for the inplace versions
M∞(v, V1, V2) = (1 + tanh((v - V1) / V2)) / 2
N∞(v, V3, V4) = (1 + tanh((v - V3) / V4)) / 2
Λ(V, V3, V4) = cosh((V - V3) / (2 * V4));
Φe(v, VSe, V0e) = 1 / (1 + exp(-VSe * (v - V0e)))
Φi(v, VSi, V0i) = 1 / (1 + exp(-VSi * (v - V0i)))
ħe(e, k_ACh) = (e^2) / (e^2 + k_ACh)
ħi(i, k_GABA) = (i^2) / (i^2 + k_GABA)

α_M(v, V7, V8, V9) = -(v - V8) / (V7 * (exp(-(v - V8) / V9) - 1))
β_M(v, V10, V11, V12) = V10 * (exp(-(v - V11) / V12))
α_H(v, V13, V14, V15) = V13 * (exp(-(v - V14) / V15))
β_H(v, V16, V17, V18) = 1 / (V16 * (exp(-(v - V17) / V18) + 1))

R∞(V, V1, V2) = (1 + tanh((V - V1) / V2)) / 2
fI(g, R, v, E) = -g * R * (v - E)

#Noise models
noise(du::Array{T2,1}, u::Array{T2,1}, p, t::T) where {T<:Real,T2} = du[end] = p[end]
noise(du::Array{T2,3}, u::Array{T2,3}, p, t::T) where {T<:Real,T2} = du[:, :, end] .= p[end]
#noise(du::CuArray{T,3}, u::CuArray{T,3}, p, t::T) where {T<:Real, T2} = du[:, :, end] .= p[end]

lansdell_noise(du, u, p, t) = du[:, :, end] .= p[end]

#ODE=============================================================================================================#
print("Loading ODE models... ")
@named SAC_Algebraic = ODESystem([
     Dt(I_ext) ~ -I_ext + I_app, #This parameter is controlled by an outside      
     Dt(v) ~ (ILeak(v) + ICa(v) + IK(v) + ITREK(v) + IACh(v) + IGABA(v) + INa(v) + I_ext + W) / C_m,
     I_Ca ~ ICa(v),
     I_Na ~ INa(v),
     I_K ~ IK(v),
     Dt(n) ~ (Λ(v) * ((N∞(v) - n))) / τn,
     Dt(m) ~ α_M(v) * (1 - m) - β_M(v) * m,
     Dt(h) ~ α_H(v) * (1 - h) - β_H(v) * h,
     Dt(c) ~ (C_0 + δ * (ICa(v)) - λ * c) / τc,
     Dt(a) ~ (α * c^4 * (1 - a) - a) / τa,
     Dt(b) ~ (β * a^4 * (1 - b) - b) / τb,
     Dt(e) ~ (ρe * Φe(v) - e) / τACh, #ρe-e, #
     Dt(i) ~ (ρi * Φi(v) - i) / τGABA, #ρi-i, #
     Dt(W) ~ -W / τw
])
ODEModel = structural_simplify(SAC_Algebraic)
println("Complete")
#SDE=============================================================================================================#
print("Loading SDE models... ")
noise_eqs = zeros(length(ODEModel.eqs))
noise_eqs[end] = 0.1

@named SDEModel = SDESystem(ODEModel, noise_eqs)
println("Complete")
#PDE=============================================================================================================#
PDEeqs = [
     Dt(v̂(x, y, t)) ~ (ÎLeak(v̂(x, y, t)) + ÎCa(v̂(x, y, t)) + ÎK(v̂(x, y, t)) + ÎTREK(v̂(x, y, t)) + ÎACh(v̂(x, y, t)) + ÎGABA(v̂(x, y, t)) + ÎNa(v̂(x, y, t)) + I_app + Ŵ(x, y, t)) / C_m,
     Dt(n̂(x, y, t)) ~ (Λ(v̂(x, y, t)) * ((N∞(v̂(x, y, t)) - n̂(x, y, t)))) / τn,
     Dt(m̂(x, y, t)) ~ α_M(v̂(x, y, t)) * (1 - m̂(x, y, t)) - β_M(v̂(x, y, t)) * m̂(x, y, t),
     Dt(ĥ(x, y, t)) ~ α_H(v̂(x, y, t)) * (1 - ĥ(x, y, t)) - β_H(v̂(x, y, t)) * ĥ(x, y, t),
     Dt(ĉ(x, y, t)) ~ (C_0 + δ * ÎCa(v̂(x, y, t)) - λ * ĉ(x, y, t)) / τc,
     Dt(â(x, y, t)) ~ (α * ĉ(x, y, t)^4 * (1 - â(x, y, t)) - â(x, y, t)) / τa,
     Dt(b̂(x, y, t)) ~ (β * â(x, y, t)^4 * (1 - b̂(x, y, t)) - b̂(x, y, t)) / τb,
     Dt(ê(x, y, t)) ~ (De * ∇²(ê(x, y, t)) + ρe*Φe(v̂(x, y, t)) - ê(x, y, t)) / τACh,
     Dt(î(x, y, t)) ~ (Di * ∇²(î(x, y, t)) + ρi*Φi(v̂(x, y, t)) - î(x, y, t)) / τGABA,
     Dt(Ŵ(x, y, t)) ~ -Ŵ(x, y, t) / τw
]