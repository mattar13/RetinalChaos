######################Modelling Toolkit for Tarchick Model
@parameters t g_leak E_leak g_Ca V1 V2 E_Ca g_K E_K g_TREK g_sAHP g_ACh k_d E_ACh I_app C_m V3 V4 τn C_0 λ δ τc α τa β τb ρ τACh k V0 σ D τw τr τs δc αs αc αr H_x
@variables v(t) n(t) c(t) a(t) b(t) e(t) W(t) r(t) s(t)
@derivatives d'~t
@register R_INF(V, VS, VH)
@register Λ(v, V3, V4)
@register Φ(v, κ, V0)
@register ħ(e, k_d)

T_model_eqs = [
          d(v)~ (
                    -g_leak*(v-E_leak)
                  + -g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
                  + -g_K*n*(v-E_K)
                  + -g_TREK*b*(v-E_K)
                  + -g_ACh*ħ(e, k_d)*(v-E_ACh)
                  + I_app
                  + W + 0*σ
                  )/C_m ,
          d(n) ~ (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn,
          d(c) ~ (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)* (v - E_Ca)) - λ*c)/τc,
          d(a) ~ (α*c^4*(1-a)-a)/τa,
          d(b) ~ (β*a^4*(1-b)-b)/τb,
          d(e) ~ (-2*D*e) + (ρ*Φ(v, k, V0) - e)/τACh,
          d(W) ~ -W/τw
]
T_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, σ]
T_ode = ODESystem(T_model_eqs)
T_sde = SDESystem(T_model_eqs, T_noise_eqs, t, T_ode.states, T_ode.ps)

#Karvourniari Model
K_model_eqs = [
          d(v) ~ (-g_leak*(v-E_leak)
                    + -g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
                    + -g_K*n*(v-E_K)
                    + -g_sAHP*r^4*(v-E_K)
                    + I_app
                    + W + 0*σ
                    )/C_m,
          d(n)~ (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn,
          d(r) ~ (αr*s*(1-r)-r)/τr,
          d(s) ~ (αs*c^4*(1-s)-s)/τs,
          d(c) ~ (-(αc/H_x)*c + C_0 + δc*(-g_Ca*R_INF(v, V1, V2)*(v-E_Ca)))/τc,
          d(W) ~ -W
]
K_noise_eqs = [0.0, 0.0, 0.0, 0.0, 0.0, σ]
K_ode = ODESystem(K_model_eqs)
K_sde = SDESystem(K_model_eqs, K_noise_eqs, t, K_ode.states, K_ode.ps)
