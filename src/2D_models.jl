#Version 0: Model without nullout
function (PDE::Network{T, :Default})(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where T
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
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p

    @. dv = (
            - g_leak*(v-E_leak)
            - g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
            - g_K*n*(v-E_K)
            - g_TREK*b*(v-E_K)
            - g_ACh*ħ(e, k_d)*(v-E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyE, PDE.My, e)
    mul!(PDE.EMx, e, PDE.Mx)
    @. PDE.DE = D*(PDE.MyE + PDE.EMx)
    @. de = PDE.DE + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 1: gTREK nullout
function (PDE::Network{T, :gHCN})(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where T
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
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p
    
    @. dv = (
            - g_leak*(v-E_leak)
            - g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
            - g_K*n*(v-E_K)
            - g_TREK*b*(v-E_K) 
            - g_ACh*ħ(e, k_d)*(v-E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyE, PDE.My, e)
    mul!(PDE.EMx, e, PDE.Mx)
    @. PDE.DE = D*(PDE.MyE + PDE.EMx)
    @. de = PDE.DE + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 2: gACh nullout
function (PDE::Network{T, :gACh})(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where T
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
            - g_leak*(v-E_leak)
            - g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
            - g_K*n*(v-E_K)
            - g_TREK*b*(v-E_K)
            - g_ACh*ħ(e, k_d)*(v-E_ACh) .* PDE.null
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyE, PDE.My, e)
    mul!(PDE.EMx, e, PDE.Mx)
    @. PDE.DE = D*(PDE.MyE + PDE.EMx)
    @. de = PDE.DE + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 3 :Acetylcholine release nullout
function (PDE::Network{T, :ρ})(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where T
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
    dW = view(dU, :,:,7)

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, D, τw, σ) = p

    @. dv = (
            - g_leak*(v-E_leak)
            - g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
            - g_K*n*(v-E_K)
            - g_TREK*b*(v-E_K)
            - g_ACh*ħ(e, k_d)*(v-E_ACh)
            + I_app
            + W
        )/C_m
    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyE, PDE.My, e)
    mul!(PDE.EMx, e, PDE.Mx)
    @. PDE.DE = D*(PDE.MyE + PDE.EMx)
    @. de = PDE.DE + (ρ*Φ(v, k, V0).*PDE.null - e)/τACh
    @. dW = -W/τw
    nothing
end

#This is the Lansdell version of the SAC model
function (PDE::Network{T, :Lansdell})(dU::AbstractArray{T}, U::AbstractArray{T}, p::AbstractArray, t::T) where T
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
            - g_leak*(v-E_leak)
            - g_Ca*R_INF(v, V1, V2)*(v-E_Ca)
            - g_K*n*(v-E_K)
            - g_TREK*b*(v-E_K)
            - g_ACh*ħ(e, k_d)*(v-E_ACh)
            + I_app
            + W
        )/C_m

    @. dr = (Λ(v, V3, V4)*(R_INF(v, V3, V4) - r) + α*s*(1-r))/τr
    @. ds = γ*Φ(v, k, V0)-s/τs
    mul!(PDE.MyE, PDE.My, a)
    mul!(PDE.EMx, a, PDE.Mx)
    @. PDE.DE = D*(PDE.MyA + PDE.AMx)
    @. da = PDE.DE + β*Φ(v, k, V0)-a/τACh
    @. dW = -W
    nothing
end
