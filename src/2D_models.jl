function ∇(du::Array{T,2}, u::Array{T, 2}, D::Float64) where T
    nx, ny = size(u)
    #These are boundary conditions for all x's at the first position
    @inbounds for y in 2:ny-1
        x = 1
        du[x,y] = D*(2u[x+1,y] + u[x,y+1] + u[x,y-1] - 4u[x,y])
    end
    
    #These are boundary conditions for all y's at the first position
    @inbounds for x in 2:nx-1
        y = 1
        du[x,y] = D*(u[x-1,y]+u[x+1,y]+ 2u[x,y+1] - 4u[x,y])
    end

    #These are boundary conditions for x's at the end
    @inbounds for y in 2:ny-1
        x = nx
        du[x,y] = D*(2u[x-1,y] + u[x,y+1] + u[x,y-1] - 4u[x,y])
    end
    
    #These are boundary conditions for all y's at the first position
    @inbounds for x in 2:ny-1
        y = ny
        du[x,y] = D*(2u[x-1,y] + u[x,y+1] + u[x,y-1] - 4u[x,y])
    end
    
    @inbounds begin
        x = 1; y = 1
        du[x,y]  = D*(2u[x+1,y]+2u[x,y+1] - 4u[x,y])
        x = 1; y = ny
        du[x,y]  = D*(2u[x+1,y]+2u[x,y-1] - 4u[x,y])
        x = nx; y = 1
        du[x,y]  = D*(2u[x-1,y]+2u[x,y+1] - 4u[x,y])
        x = nx; y = ny
        du[x,y]  = D*(2u[x-1,y]+2u[x,y-1] - 4u[x,y])
    end
    
    @inbounds for x in 2:nx-1, y in 2:ny-1
        du[x,y] = D*(u[x-1,y]+u[x+1,y]+u[x,y-1]+u[x,y+1]-4u[x,y])
    end
    du
end
#Version -1: Testing Unwinding and Inbounds
function (PDE::Network{Array{T,2}, :Unwinding})(dU::Array{T,3}, U::Array{T,3}, p::Array{T,1}, t::T) where T <: Real
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

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              I_n(g_leak,  1.0, v, E_leak)
            + I_n(g_Ca, R_INF(v, V1, V2), v, E_Ca)
            + I_n(g_K, n , v, E_K)
            + I_n(g_TREK, b, v, E_K)
            + I_n(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    #mul!(PDE.MyA, PDE.My, e)
    #mul!(PDE.AMx, e, PDE.Mx)
    #PDE.DA = ∇(PDE.DA, e, D)
    @. de = ∇(de, e, D) + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 0: Model without nullout
function (PDE::Network{Array{T,2}, :Default})(dU::Array{T,3}, U::Array{T,3}, p::Array{T,1}, t::T) where T <: Real
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

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              I_n(g_leak,  1.0, v, E_leak)
            + I_n(g_Ca, R_INF(v, V1, V2), v, E_Ca)
            + I_n(g_K, n , v, E_K)
            + I_n(g_TREK, b, v, E_K)
            + I_n(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 1: gTREK nullout
function (PDE::Network{Array{T,2}, :gTREK})(dU::Array{T,3}, U::Array{T,3}, p::Array{T,1}, t::T) where T <: Real
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

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              I_n(g_leak,  1.0, v, E_leak)
            + I_n(g_Ca, R_INF(v, V1, V2), v, E_Ca)
            + I_n(g_K, n , v, E_K)
            + I_n(g_TREK, b, v, E_K) .* PDE.null
            + I_n(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 2: gACh nullout
function (PDE::Network{Array{T,2}, :gACh})(dU::Array{T,3}, U::Array{T,3}, p::Array{T,1}, t::T) where T <: Real
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

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              I_n(g_leak,  1.0, v, E_leak)
            + I_n(g_Ca, R_INF(v, V1, V2), v, E_Ca)
            + I_n(g_K, n , v, E_K)
            + I_n(g_TREK, b, v, E_K)
            + I_n(g_ACh, ħ(e, k_d), v, E_ACh) .* PDE.null
            + I_app
            + W
        )/C_m

    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0) - e)/τACh
    @. dW = -W/τw
    nothing
end

#Version 3 :Acetylcholine release nullout
function (PDE::Network{Array{T,2}, :ρ})(dU::Array{T,3}, U::Array{T,3}, p::Array{T,1}, t::T) where T <: Real
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

    (g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_d, E_ACh, I_app, C_m, V3, V4, τn, C_0, λ, δ, τc, α, τa, β, τb, ρ, τACh, k, V0, σ, D, τw) = p

    @. dv = (
              I_n(g_leak,  1.0, v, E_leak)
            + I_n(g_Ca, R_INF(v, V1, V2), v, E_Ca)
            + I_n(g_K, n , v, E_K)
            + I_n(g_TREK, b, v, E_K)
            + I_n(g_ACh, ħ(e, k_d), v, E_ACh)
            + I_app
            + W
        )/C_m
    @. dn = (Λ(v, V3, V4) * ((R_INF(v, V3, V4) - n)))/τn
    @. dc = (C_0 + δ*(-g_Ca*R_INF(v, V1, V2)*(v - E_Ca)) - λ*c)/τc
    @. da =  (α*c^4*(1-a) - a)/τa
    @. db =  (β*a^4*(1-b) - b)/τb
    mul!(PDE.MyA, PDE.My, e)
    mul!(PDE.AMx, e, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. de = PDE.DA + (ρ*Φ(v, k, V0)*PDE.null - e)/τACh
    @. dW = -W/τw
    nothing
end

lansdell_pars = [:I_app, :E_Ca, :E_K, :E_Leak, :E_ACh, :V1, :V2, :V3, :V4, :g_Ca, :g_K, :g_Leak, :λ, :g_ACh, :δ, :C_m, :τr, :τs, :τACh, :γ, :α, :β, :k, :V0, :D, :μ]
lansdell_conds = [:v, :r, :s, :a, :W]

#This is the Lansdell version of the SAC model
function (PDE::Network{Array{T,2}, :Lansdell})(dU::Array{T,3}, U::Array{T,3}, p::Array{T,1}, t::T) where T <: Real
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
              I_n(g_leak,  1.0, v, E_leak)
            + I_n(g_Ca, R_INF(v, V1, V2), v, E_Ca)
            + I_n(g_K, r , v, E_K)
            + I_n(g_ACh, h(a, δ), v, E_ACh)
            + fIn(abs(W), λ, v, E_Ca)
            + I_app
        )/C_m

    @. dr = (Λ(v, V3, V4)*(R_INF(v, V3, V4) - r) + α*s*(1-r))/τr
    @. ds = γ*Φ(v, k, V0)-s/τs
    mul!(PDE.MyA, PDE.My, a)
    mul!(PDE.AMx, a, PDE.Mx)
    @. PDE.DA = D*(PDE.MyA + PDE.AMx)
    @. da = PDE.DA + β*Φ(v, k, V0)-a/τACh
    @. dW = -W
    nothing
end
