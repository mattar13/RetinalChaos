function ∇(du, u, D)
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
        du[x,y] = D*(u[x-1,y] + u[x+1,y] + 2u[x,y-1] - 4u[x,y])
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
end

mutable struct Network{T, N} <: Function
      Mx::AbstractArray{T,2}
      My::AbstractArray{T,2}
      MyE::AbstractArray{T,2}
      EMx::AbstractArray{T,2}
      DE::AbstractArray{T,2}
      null::AbstractArray{T,2}
end

"""
This is a constructor for the Network object with a version flag option
"""
Network(Mx::AbstractArray{T}, My::AbstractArray{T}, MyE::AbstractArray{T}, EMx::AbstractArray{T}, DE::AbstractArray{T}, null::AbstractArray{T}, null_param::Symbol) where T  = Network{T, null_param}(Mx, My, MyE, EMx, DE, null)

"""
This constructs the PDE function so that it can be called
"""
function Network(nx::Real, ny::Real; 
            gpu::Bool = false, μ::T = 0.75, version = :gACh,
            DX::Tuple{Float64, Float64} = (-2.0, 1.0), DY::Tuple{Float64, Float64} = (-2.0, 1.0)
        )where T <: Real
    #We want to ensure both nx and ny are integers    
    nx = nx |> Int64 
    ny = ny |> Int64
    x_dv = repeat([DX[1]], nx)
    x_uv = repeat([DX[2]], nx-2)
    x_lv = repeat([DX[2]], nx-2)
    push!(x_uv, abs(DX[1]))
    pushfirst!(x_lv, abs(DX[1]))

    #Set up y diffusion steps
    y_dv = repeat([DY[1]], ny)
    y_uv = repeat([DY[2]], ny-2)
    y_lv = repeat([DY[2]], ny-2)
    push!(y_lv, abs(DY[1]))
    pushfirst!(y_uv, abs(DY[1]))
    Mx = Array(Tridiagonal(x_lv, x_dv, x_uv))
    My = Array(Tridiagonal(y_lv, y_dv, y_uv))
    MyE = zeros(ny, nx);
    EMx = zeros(ny, nx);
    DE = zeros(ny, nx);
    d = Binomial(1, μ|>Float64)
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

function diffuse(lattice, D, PDE::Network)
    mul!(PDE.MyE, PDE.My, lattice)
    mul!(PDE.EMx, lattice, PDE.Mx)
    DE = D*(PDE.MyE + PDE.EMx)
    return lattice + DE
end
