mutable struct Network{A, N} <: Function
      Mx::Union{Tridiagonal, A}
      My::Union{Tridiagonal, A}
      MyA::A
      AMx::A
      DA::A
      null::A
end

"""
This is a constructor for the Network object with a version flag option
"""
Network(Mx, My, MyA, AMx, DA, null, null_param::Symbol) = Network{typeof(MyA), null_param}(Mx, My, MyA, AMx, DA, null)

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

"""
This constructs the PDE function so that it can be called
"""
function Network(nx::Int64, ny::Int64; gpu::Bool = false, μ::Float64 = 0.75, version = :gACh,
        DX::Tuple{Float64, Float64} = (-2.0, 1.0), DY::Tuple{Float64, Float64} = (-2.0, 1.0))
    #Set up x diffusion steps
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
    Mx = Tridiagonal(x_lv, x_dv, x_uv)
    My = Tridiagonal(y_lv, y_dv, y_uv)
    MyA = zeros(ny, nx);
    AMx = zeros(ny, nx);
    DA = zeros(ny, nx);
    d = Binomial(1, μ)
    null = Float64.(rand(d, (ny, nx)))
    if gpu
        gMx = CuArray(Float32.(Mx))
        gMy = CuArray(Float32.(Mx))
        gAMx = CuArray(Float32.(Mx))
        gMyA = CuArray(Float32.(Mx))
        gDA = CuArray(Float32.(Mx))
        gNull = CuArray(Float32.(null))
        return Network(gMx, gMy, gMyA, gAMx, gDA, gNull, version)
    else
        return Network(Mx, My, MyA, AMx, DA, null, version)
    end
end

function diffuse(lattice, D, PDE::Network)
    mul!(PDE.MyA, PDE.My, lattice)
    mul!(PDE.AMx, lattice, PDE.Mx)
    DA = D*(PDE.MyA + PDE.AMx)
    return lattice + DA
end
