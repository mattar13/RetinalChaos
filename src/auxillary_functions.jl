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

∇²(u) = Dxx(u) + Dyy(u) #This is the diffusion aspect of the equation
