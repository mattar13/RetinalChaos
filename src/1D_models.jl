function OLD_ODE(dU::AbstractArray{T2}, U::AbstractArray{T2}, p::AbstractArray{T}, t::T) where {T<:Real,T2}
     #Extract all of the parameters first
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

     @. dv = (
          -g_leak * (v - E_leak)
          -
          g_Ca * R_INF(v, V1, V2) * (v - E_Ca)
          -
          g_K * n * (v - E_K)
          -
          g_TREK * b * (v - E_K)
          -
          g_ACh * ħ(e, k_d) * (v - E_ACh)
          +
          I_app
     ) / C_m

     @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
     @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
     @. da = (α * c^4 * (1 - a) - a) / τa
     @. db = (β * a^4 * (1 - b) - b) / τb
     #@. de = (ρ*Φ(v, k, V0) - e)/τACh
     @. de = (ρe - e) #Keep release constant
     @. dW = -W / τw
     #nothing
     dU
end

function OLD_SDE(dU::AbstractArray{T2}, U::AbstractArray{T2}, p::AbstractArray{T}, t::T) where {T<:Real,T2}
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

     @.dv = (
          -g_leak * (v - E_leak)
          -
          g_Ca * R_INF(v, V1, V2) * (v - E_Ca)
          -
          g_K * n * (v - E_K)
          -
          g_TREK * b * (v - E_K)
          -
          g_ACh * ħ(e, k_d) * (v - E_ACh)
          + I_app
          + W
     ) / C_m

     @.dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
     @.dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
     @.da = (α * c^4 * (1 - a) - a) / τa
     @.db = (β * a^4 * (1 - b) - b) / τb
     #@.de = (ρ*Φ(v, k, V0) - e)/τACh
     @.de = (ρe - e) #Keep release constant
     @.dW = -W / τw
     #nothing
     dU
end

#This includes the GABA term
function T_ODE(dU::AbstractArray{T2}, U::AbstractArray{T2}, p::AbstractArray{T}, t::T) where {T<:Real,T2}
     #Extract all of the parameters first
     (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_ACh, E_ACh, g_GABA, k_GABA, E_GABA, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρe, ρi, τACh, τGABA, ke, ki, V0e, V0i, De, Di, τw, σ) = p
     v = view(U, 1)
     n = view(U, 2)
     c = view(U, 3)
     a = view(U, 4)
     b = view(U, 5)
     e = view(U, 6)
     i = view(U, 7)
     W = view(U, 8)

     dv = view(dU, 1)
     dn = view(dU, 2)
     dc = view(dU, 3)
     da = view(dU, 4)
     db = view(dU, 5)
     de = view(dU, 6)
     di = view(dU, 7)
     dW = view(dU, 8)

     @. dv = (-g_leak * (v - E_leak) - g_Ca * R_INF(v, V1, V2) * (v - E_Ca) - g_K * n * (v - E_K) - g_TREK * b * (v - E_K) - g_ACh * ħ(e, k_ACh) * (v - E_ACh) - g_GABA * ħ(i, k_GABA) * (v - E_GABA) + I_app) / C_m
     @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
     @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
     @. da = (α * c^4 * (1 - a) - a) / τa
     @. db = (β * a^4 * (1 - b) - b) / τb
     #@. de = (ρe * Φ(v, ke, V0e) - e) / τACh
     #@. di = (ρi * Φ(v, ki, V0i) - i) / τGABA
     @. de = (ρe-e) #Keep release constant
     @. di = (ρi-i) #Keep release constant
     @. dW = -W / τw
     #nothing
     dU
end

#This includes the GABA term
function T_SDE(dU::AbstractArray{T2}, U::AbstractArray{T2}, p::AbstractArray{T}, t::T) where {T<:Real,T2}
     #Extract all of the parameters first
     (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_ACh, E_ACh, g_GABA, k_GABA, E_GABA, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρe, ρi, τACh, τGABA, ke, ki, V0e, V0i, De, Di, τw, σ) = p
     v = view(U, 1)
     n = view(U, 2)
     c = view(U, 3)
     a = view(U, 4)
     b = view(U, 5)
     e = view(U, 6)
     i = view(U, 7)
     W = view(U, 8)

     dv = view(dU, 1)
     dn = view(dU, 2)
     dc = view(dU, 3)
     da = view(dU, 4)
     db = view(dU, 5)
     de = view(dU, 6)
     di = view(dU, 7)
     dW = view(dU, 8)

     @. dv = (-g_leak * (v - E_leak) - g_Ca * R_INF(v, V1, V2) * (v - E_Ca) - g_K * n * (v - E_K) - g_TREK * b * (v - E_K) - g_ACh * ħ(e, k_ACh) * (v - E_ACh) - g_GABA * ħ(i, k_GABA) * (v - E_GABA) + I_app + W) / C_m
     @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
     @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
     @. da = (α * c^4 * (1 - a) - a) / τa
     @. db = (β * a^4 * (1 - b) - b) / τb
     #@. de = (ρe * Φ(v, ke, V0e) - e) / τACh
     #@. di = (ρi * Φ(v, ki, V0i) - i) / τGABA
     #@. de = ρe*abs(sin(t/20000)) - e#(ρe-e) #Sine wave
     #@. di = ρi*abs(sin(t/20000)) - i #(ρi-i) #Sine wave
     @. de = (ρe-e) #Keep release constant
     @. di = (ρi-i) #Keep release constant
     @. dW = -W / τw
     #nothing
     dU
end