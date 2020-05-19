######################Modelling Toolkit for Tarchick Model
@parameters t g_leak E_leak g_Ca V1 V2 E_Ca g_K E_K g_TREK g_ACh k_d E_ACh I_app C_m V3 V4 τn C_0 λ δ τc α τa β τb ρ τACh k V0 σ D τw
all_pars = [g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw]
@variables v(t) n(t) c(t) a(t) b(t) ACh(t) W(t)
@derivatives d'~t

T_model_eqs = [d(v)~ (
          -g_leak*(v-E_leak)
        + -g_Ca*M_INF(v, V1, V2)*(v-E_Ca)
        + -g_K*n*(v-E_K)
        + -g_TREK*b*(v-E_K)
        + -g_ACh*ħ(ACh, k_d)*(v-E_ACh)
        + I_app
        + W
        )/C_m ,
        d(n) ~ (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn,
        d(c) ~ (C_0 + δ*(-g_Ca*M_INF(v, V1, V2)* (v - E_Ca)) - λ*c)/τc,
        d(a) ~ (α*c^4*(1-a)-a)/τa,
        d(b) ~ (β*a^4*(1-b)-b)/τb,
        d(ACh) ~ (-2*D*ACh) + (ρ*Φ(v, k, V0) - ACh)/τACh,
        d(W) ~ -W/τw
]
T_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, σ]
T_model_sys = SDESystem(T_model_eqs, T_noise_eqs, t, [v, n, c, a, b, ACh, W], all_pars)

#Karvourniari Model
K_pars = [I_app,C_m,g_leak,g_Ca,g_K,g_sAHP,E_leak,E_Ca,E_K,V1,V2,V3,V4,τn,τr,τs,τc,δc,αs,αc,αr,H_x,C_0,σ]
K_model_eqs = [
  d(v) ~ (-g_leak*(v-E_leak)
        + -g_Ca*M_INF(v, V1, V2)*(v-E_Ca)
        + -g_K*n*(v-E_K)
        + -g_TREK*r^4*(v-E_K)
        + -g_ACh*ħ(ACh, k_d)*(v-E_ACh)
        + I_app
        + W
        )/C_m,
  d(n) ~ (Λ(v, V3, V4) * ((N_INF(v, V3, V4) - n)))/τn,
  d(r) ~ (αr*s*(1-r)-r)/τr,
  d(s) ~ (αs*c^4*(1-s)-s)/τs,
  d(c) ~ -(αc/H_x)*c + C_0 + δc*(-g_Ca*M_INF(v, V1, V2)*(v-E_Ca)),
  d(W) ~ -W
]
K_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, σ]
K_sys = SDESystem(K_model_eqs, K_noise_eqs, t, [v, n, r, s, c, W], K_pars)
