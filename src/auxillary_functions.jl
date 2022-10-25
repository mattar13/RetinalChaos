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