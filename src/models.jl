#Auxillary functions
M_INF(V, V1, V2) = (1 + tanh((V - V1)/V2))/2;
N_INF(V, V3, V4) = (1 + tanh((V - V3)/V4))/2;
Λ(V, V3, V4) = cosh((V-V3)/(2*V4));
Φ(v, κ, V_0) = 1/(1 + exp(-κ * (v - V_0)))
ħ(a, K_d) = (a^2)/(a^2 + K_d)
h(a, δ) = (δ*a^2)/(1+δ*a^2)

CuArrays.@cufunc M_INF(V, V1, V2) = (1 + tanh((V - V1)/V2))/2;
CuArrays.@cufunc N_INF(V, V3, V4) = (1 + tanh((V - V3)/V4))/2;
CuArrays.@cufunc Λ(V, V3, V4) = cosh((V-V3)/(2*V4));
CuArrays.@cufunc Φ(v, κ, V_0) = 1/(1 + exp(-κ * (v - V_0)))
CuArrays.@cufunc ħ(a, K_d) = (a^2)/(a^2 + K_d)
CuArrays.@cufunc h(a, δ) = (δ*a^2)/(1+δ*a^2)

M_INF(V::CuArray, V1, V2) = M_INF.(V, V1, V2)
N_INF(V::CuArray, V3, V4) = N_INF.(V, V3, V4)
Λ(V::CuArray, V3, V4) = LAM.(V, V3, V4) 
Φ(v::CuArray, κ, V_0) = Φ.(v, κ, V_0)
ħ(a::CuArray, K_d) = ħ.(a, K_d)
h(a::CuArray, δ) = h.(a, δ)

#This is the generalized current calculation function
fI(g::Float64, r, v, e::Float64) = -g*r*(v-e)
#This is for calculating filtered shot noise
fIn(W::Float64, λ::Float64, v, e::Float64) = -(-log(W)/λ)*(v-e)  


#p_find(p) = findall(x -> x == p, BurstModel.params)
#u_find(u) = findall(x -> x == u, BurstModel.syms)

include("pde_models.jl")
include("1D_models.jl") 
include("2D_models.jl")    



#Noise models
noise(du, u, p, t) = du[end] = p[30]
noise_2D(du, u, p, t) = du[:,:,end] .= p[30]
lansdell_noise(du, u, p, t) = du[:,:,end] .= p[end]

