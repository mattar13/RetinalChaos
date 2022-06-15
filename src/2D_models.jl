#Version 0: Model without nullout
function T_PDE(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where {T}
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :, :, 7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, Vs, V0, D, τw, σ) = p

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
        + I_app
        + W
    ) / C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
    @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
    @. da = (α * c^4 * (1 - a) - a) / τa
    @. db = (β * a^4 * (1 - b) - b) / τb
    @. de = (ρ * Φ(v, k, V0) - e) / τACh
    ∇(de, e, D; dX=dX, dY=dY) #This is the diffusion step
    @. dW = -W / τw
    nothing
end

#Version 1: gTREK nullout
function gHCN_PDE(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where {T}
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :, :, 7)

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
        + I_app
        + W
    ) / C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
    @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
    @. da = (α * c^4 * (1 - a) - a) / τa
    @. db = (β * a^4 * (1 - b) - b) / τb

    @. de = (ρ * Φ(v, k, V0) - e) / τACh
    ∇(de, e, D) #This is the diffusion step
    @. dW = -W / τw
    nothing
end

#Version 2: gACh nullout
function gACh_PDE(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T, null::AbstractArray{T}) where {T}
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :, :, 7)

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
        g_ACh * ħ(e, k_d) * (v - E_ACh) .* null
        + I_app
        + W
    ) / C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
    @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
    @. da = (α * c^4 * (1 - a) - a) / τa
    @. db = (β * a^4 * (1 - b) - b) / τb
    @. de = (ρ * Φ(v, k, V0) - e) / τACh
    ∇(de, e, D) #This is the diffusion step
    @. dW = -W / τw
    nothing
end

#Version 3 :Acetylcholine release nullout
function Rho_PDE(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where {T}
    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    W = view(U, :, :, 7)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    dW = view(dU, :, :, 7)

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
        + I_app
        + W
    ) / C_m
    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
    @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
    @. da = (α * c^4 * (1 - a) - a) / τa
    @. db = (β * a^4 * (1 - b) - b) / τb
    @. de = (ρ * Φ(v, k, V0) .* null - e) / τACh
    ∇(de, e, D) #This is the diffusion step
    @. dW = -W / τw
    nothing
end

#Version 4: GABA sensitive model.
function GABA_PDE(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T; dXe=(1.0, 1.0), dYe=(1.0, 1.0), dXi = (1.0, 1.0), dYi = (0.1, 1.9)) where {T}
    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_ACh, E_ACh, g_GABA, k_GABA, E_GABA, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρe, ρi, τACh, τGABA, Vse, Vsi, V0e, V0i, De, Di, τw, σ) = p

    v = view(U, :, :, 1)
    n = view(U, :, :, 2)
    c = view(U, :, :, 3)
    a = view(U, :, :, 4)
    b = view(U, :, :, 5)
    e = view(U, :, :, 6)
    i = view(U, :, :, 7)
    W = view(U, :, :, 8)

    dv = view(dU, :, :, 1)
    dn = view(dU, :, :, 2)
    dc = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    db = view(dU, :, :, 5)
    de = view(dU, :, :, 6)
    di = view(dU, :, :, 7)
    dW = view(dU, :, :, 8)

    @. dv = (-g_leak * (v - E_leak)
             -g_Ca * R_INF(v, V1, V2) * (v - E_Ca)
             -g_K * n * (v - E_K)
             -g_TREK * b * (v - E_K)
             -g_ACh * ħ(e, k_ACh) * (v - E_ACh)
             -g_GABA * ħ(i, k_GABA) * (v - E_GABA)
             + I_app
             + W
    ) / C_m
    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n))) / τn
    @. dc = (C_0 + δ * (-g_Ca * R_INF(v, V1, V2) * (v - E_Ca)) - λ * c) / τc
    @. da = (α * c^4 * (1 - a) - a) / τa
    @. db = (β * a^4 * (1 - b) - b) / τb

    @. de = (ρe * Φ(v, Vse, V0e) - e) / τACh
    ∇(de, e, De; dX = dXe, dY = dYe) #This is the diffusion step NOTE there is an error with the naming of x and y i need to fix

    @. di = (ρi * Φ(v, Vsi, V0i) - i) / τGABA
    ∇(di, i, Di; dX = dXi, dY = dYi) #This is the diffusion step

    @. dW = -W / τw
    nothing
end

#This is the Lansdell version of the SAC model
function Lansdell_PDE(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where {T}
    v = view(U, :, :, 1)
    r = view(U, :, :, 2)
    s = view(U, :, :, 3)
    a = view(U, :, :, 4)
    W = view(U, :, :, 5)

    dv = view(dU, :, :, 1)
    dr = view(dU, :, :, 2)
    ds = view(dU, :, :, 3)
    da = view(dU, :, :, 4)
    dW = view(dU, :, :, 5)

    (I_app, E_Ca, E_K, E_leak, E_ACh, V1, V2, V3, V4, g_Ca, g_K, g_leak, λ, g_ACh, δ, C_m, τr, τs, τACh, γ, α, β, k, V0, D, μ) = p

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
        + I_app
        + W
    ) / C_m

    @. dr = (Λ(v, V3, V4) * (R_INF(v, V3, V4) - r) + α * s * (1 - r)) / τr
    @. ds = γ * Φ(v, k, V0) - s / τs
    mul!(PDE.MyE, PDE.My, a)
    mul!(PDE.EMx, a, PDE.Mx)
    @. PDE.DE = D * (PDE.MyA + PDE.AMx)
    @. da = PDE.DE + β * Φ(v, k, V0) - a / τACh
    @. dW = -W
    nothing
end
