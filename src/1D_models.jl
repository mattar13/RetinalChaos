######################Modelling Toolkit for Tarchick Model
@parameters t g_leak E_leak g_Ca V1 V2 E_Ca g_K E_K g_TREK g_ACh k_d E_ACh I_app C_m V3 V4 τn C_0 λ δ τc α τa β τb ρ τACh k V0 σ D τw
all_pars = [g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw]
@variables v(t) n(t) c(t) a(t) b(t) ACh(t) W(t)
@derivatives d'~t

TM_eqs = [d(v)~ (
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
noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, σ]
TModel_sys = SDESystem(TM_eqs, noise_eqs, t, [v, n, c, a, b, ACh, W], all_pars)
