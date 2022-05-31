"""
This stencil function saves memory
"""

function ∇(du, u, D; dX = (1.0, 1.0), dY = (1.0, 1.0))
    nx, ny = size(u)
    #These are boundary conditions for all x's at the first position
    sum_dx = sum(dX)
    sum_dy = sum(dY)
    sum_d = sum_dx + sum_dy
    @inbounds for y in 2:ny-1
        x = 1
        du[x, y] += D * (
            sum_dx * u[x+1, y] + 
            dY[1] * u[x, y+1] + 
            dY[2] * u[x, y-1] - 
            sum_d * u[x, y])
    end

    #These are boundary conditions for all y's at the first position
    @inbounds for x in 2:nx-1
        y = 1
        du[x, y] += D * (
            dX[1] * u[x+1, y] + 
            dX[2] * u[x-1, y] + 
            sum_dy * u[x, y+1] - 
            sum_d * u[x, y])
    end

    #These are boundary conditions for x's at the end
    @inbounds for y in 2:ny-1
        x = nx
        du[x, y] += D * (
            sum_dx * u[x-1, y] + 
            dY[1] * u[x, y+1] + 
            dY[2] * u[x, y-1] - 
            sum_d * u[x, y])
    end

    #These are boundary conditions for all y's at end position
    @inbounds for x in 2:ny-1
        y = ny
        du[x, y] += D * (
            dX[1] * u[x+1, y] + 
            dX[2] * u[x-1, y] + 
            sum_dy * u[x, y-1] - 
            sum_d * u[x, y])
    end

    @inbounds begin
        x = 1
        y = 1
        du[x, y] += D * (sum_dx * u[x+1, y] + sum_dy * u[x, y+1] - sum_d * u[x, y])
        x = 1
        y = ny
        du[x, y] += D * (sum_dx * u[x+1, y] + sum_dy * u[x, y-1] - sum_d * u[x, y])
        x = nx
        y = 1
        du[x, y] += D * (sum_dx * u[x-1, y] + sum_dy * u[x, y+1] - sum_d * u[x, y])
        x = nx
        y = ny
        du[x, y] += D * (sum_dx * u[x-1, y] + sum_dy * u[x, y-1] - sum_d * u[x, y])
    end

    @inbounds for x in 2:nx-1, y in 2:ny-1
        du[x, y] += D * (
            dX[1] * u[x+1, y] + 
            dX[2] * u[x-1, y] + 
            dY[1] * u[x, y-1] + 
            dY[2] * u[x, y+1] - 
            sum_d * u[x, y])
    end
end

#= #We may switch to something that takes less memory like an unrolled stencil
mutable struct Network{T,N} <: Function
    Mx::AbstractArray{T,2}
    My::AbstractArray{T,2}
    MyE::AbstractArray{T,2}
    EMx::AbstractArray{T,2}
    DE::AbstractArray{T,2}
    DI::AbstractArray{T,2}
    null::AbstractArray{T,2}
end


struct Stencil{T}
    Mx::Tridiagonal{T,Vector{T}}
    My::Tridiagonal{T,Vector{T}}
    AMx::Array{T}
    MyA::Array{T}
end


function diffusion_stencil(nx::Int64, ny::Int64; dX::Tuple{Float64,Float64} = (1.0, 1.0), dY::Tuple{Float64,Float64} = (1.0, 1.0)) where {T}
    x_CENTER = repeat([-sum(dX)], nx)
    x_RIGHT = repeat([dX[1]], nx - 2)
    x_LEFT = repeat([dX[2]], nx - 2)
    pushfirst!(x_RIGHT, abs(sum(dX)))
    push!(x_LEFT, abs(sum(dX)))

    y_CENTER = repeat([-sum(dY)], ny)
    y_UP = repeat([dY[1]], ny - 2)
    y_DOWN = repeat([dY[2]], ny - 2)
    pushfirst!(y_UP, abs(sum(dY)))
    push!(y_DOWN, abs(sum(dY)))

    Mx = Tridiagonal(x_LEFT, x_CENTER, x_RIGHT) #Keeping these as stencils helps keep memory usage low
    My = Tridiagonal(y_UP, y_CENTER, y_DOWN)
    return Stencil(Mx, My, zeros(nx, ny), zeros(ny, nx))
end

function diffuse(net::Stencil{T}, A::Array{T}, D::T) where T
    mul!(net.MyA, net.My, A)
    mul!(net.AMx, A, net.Mx)
    DA = D * (net.MyA + net.AMx)
    return DA + A
end



"""
This is a constructor for the Network object with a version flag option
"""
Network(Mx::AbstractArray{T}, My::AbstractArray{T}, MyE::AbstractArray{T}, EMx::AbstractArray{T}, DE::AbstractArray{T}, null::AbstractArray{T}, null_param::Symbol) where {T} = Network{T,null_param}(Mx, My, MyE, EMx, DE, null)


"""
This constructs the PDE function so that it can be called
"""
function Network(nx::Real, ny::Real;
    gpu::Bool = false, μ::T = 0.75, version = :gACh,
    DX::Tuple{Float64,Float64} = (-2.0, 1.0), DY::Tuple{Float64,Float64} = (-2.0, 1.0)
) where {T<:Real}
    #We want to ensure both nx and ny are integers    
    nx = nx |> Int64
    ny = ny |> Int64
    x_dv = repeat([DX[1]], nx)
    x_uv = repeat([DX[2]], nx - 2)
    x_lv = repeat([DX[2]], nx - 2)
    push!(x_uv, abs(DX[1]))
    pushfirst!(x_lv, abs(DX[1]))

    #Set up y diffusion steps
    y_dv = repeat([DY[1]], ny)
    y_uv = repeat([DY[2]], ny - 2)
    y_lv = repeat([DY[2]], ny - 2)
    push!(y_lv, abs(DY[1]))
    pushfirst!(y_uv, abs(DY[1]))
    Mx = Array(Tridiagonal(x_lv, x_dv, x_uv))
    My = Array(Tridiagonal(y_lv, y_dv, y_uv))
    MyE = zeros(ny, nx)
    EMx = zeros(ny, nx)
    DE = zeros(ny, nx)
    d = Binomial(1, μ |> Float64)
    null = Float64.(rand(d, (ny, nx)))
    if gpu
        #Float32 is the most efficient form of CUDA arrays
        gMx = CuArray(Float32.(Mx))
        gMy = CuArray(Float32.(Mx))
        gEMx = CuArray(Float32.(Mx))
        gMyE = CuArray(Float32.(Mx))
        gDE = CuArray(Float32.(Mx))
        gNull = CuArray(Float32.(null))
        return Network(gMx, gMy, gMyE, gEMx, gDE, gNull, version)
    else
        return Network(Mx, My, MyE, EMx, DE, null, version)
    end
end

=#