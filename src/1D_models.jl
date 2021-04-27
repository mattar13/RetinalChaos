function T_ode(dU::AbstractArray{T2}, U::AbstractArray{T2}, p::AbstractArray{T}, t::T) where {T <: Real, T2}
	#Extract all of the parameters first
	(g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p
	v = view(U, 1)
	n = view(U, 2)
	c = view(U, 3)
	a = view(U, 4)
	b = view(U, 5)
	e = view(U, 6)
	W = view(U, 7)

	dv = view(dU, 1)
	dn = view(dU, 2)
	dc = view(dU, 3)
	da = view(dU, 4)
	db = view(dU, 5)
	de = view(dU, 6)
	dW = view(dU, 7)

	(g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p

	dv = (
			- g_leak*(v-E_leak)
			- g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
			- g_K*n*(v-E_K)
			- g_TREK*b*(v-E_K)
			- g_ACh*ħ(e, k_d)*(v-E_ACh)
			+ I_app
	)/C_m

	dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
	dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
	da = (α*c^4*(1-a) - a)/τa
	db = (β*a^4*(1-b) - b)/τb
	de = (ρ*Φ(v, k, V0) - e)/τACh
	dW = -W/τw
	#nothing
	dU
end

function T_sde(dU::AbstractArray{T2}, U::AbstractArray{T2}, p::AbstractArray{T}, t::T) where {T <: Real, T2}
	#Extract all of the parameters first
	(g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p
	v = view(U, 1)
	n = view(U, 2)
	c = view(U, 3)
	a = view(U, 4)
	b = view(U, 5)
	e = view(U, 6)
	W = view(U, 7)

	dv = view(dU, 1)
	dn = view(dU, 2)
	dc = view(dU, 3)
	da = view(dU, 4)
	db = view(dU, 5)
	de = view(dU, 6)
	dW = view(dU, 7)

	(g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p

	dv = (
			- g_leak*(v-E_leak)
			- g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
			- g_K*n*(v-E_K)
			- g_TREK*b*(v-E_K)
			- g_ACh*ħ(e, k_d)*(v-E_ACh)
			+ I_app
			+ W
	)/C_m

	dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
	dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
	da = (α*c^4*(1-a) - a)/τa
	db = (β*a^4*(1-b) - b)/τb
	de = (ρ*Φ(v, k, V0) - e)/τACh
	dW = -W/τw
	#nothing
	dU
end



#Modellling toolkit may not work so well
######################Modelling Toolkit for Tarchick Model
#= #Register parameters for ionic currents 
@parameters g_leak E_leak #Leaky currents
@parameters g_Ca V1 V2 E_Ca #Cav currents
@parameters g_K V3 V4 E_K #Kv currents
@parameters g_TREK g_sAHP #TREK and sAHP currents
@parameters g_ACh k_d E_ACh #Acetylcholine current parameters
@parameters I_app C_m  #Extras
@parameters τn τc τa τb τACh τw τr τs #Time constants
@parameters C_0 V0
@parameters λ δ α β ρ k σ D #Relationships in my model
@parameters δc αs αc αr H_x VS VH g_n R E_n #Relationships in other models
#Register variables
@variables t v(t) n(t) c(t) a(t) b(t) e(t) W(t) r(t) s(t)
#@derivatives d'~ Differential(t) #This is deprecated
d = Differential(t)

T_model_eqs = [
					d(v)~ (
										-g_leak*(v-E_leak)
									+ -g_Ca*M_INF(v, V1, V2)*(v-E_Ca)
									+ -g_K*n*(v-E_K)
									+ -g_TREK*b*(v-E_K)
									+ -g_ACh*ħ(e, k_d)*(v-E_ACh)
									+ I_app
									+ W + 0*σ #This 
									)/C_m ,
					d(n) ~ (Λ(v, V3, V4) * (N_INF(v, V3, V4) - n))/τn,
					d(c) ~ (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc,
					d(a) ~ (α*c^4*(1-a)-a)/τa,
					d(b) ~ (β*a^4*(1-b)-b)/τb,
					d(e) ~ (-2*D*e) + (ρ*Φ(v, k, V0) - e)/τACh,
					d(W) ~ -W/τw
]
T_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, σ]
T_ode = ODESystem(T_model_eqs)
T_sde = SDESystem(T_model_eqs, T_noise_eqs, t, T_ode.states, T_ode.ps)
sym_ps = map(x -> x |> Symbol, T_ode.ps)
sym_cs = map(x -> x |> Symbol, T_ode.states)

#Karvourniari Model
K_model_eqs = [
					d(v) ~ (-g_leak*(v-E_leak)
										+ -g_Ca*M_INF(v, V1, V2)*(v-E_Ca)
										+ -g_K*n*(v-E_K)
										+ -g_sAHP*r^4*(v-E_K)
										+ I_app
										+ W + 0*σ
										)/C_m,
					d(n)~ (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn,
					d(r) ~ (αr*s*(1-r)-r)/τr,
					d(s) ~ (αs*c^4*(1-s)-s)/τs,
					d(c) ~ (-(αc/H_x)*c + C_0 + δc*(-g_Ca*M_INF(v, V1, V2)*(v-E_Ca)))/τc,
					d(W) ~ -W
]
K_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, σ]
K_ode = ODESystem(K_model_eqs)
K_sde = SDESystem(K_model_eqs, K_noise_eqs, t, K_ode.states, K_ode.ps)
 =#