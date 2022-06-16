"""
This stencil function saves memory
"""

function ∇(du, u, D; dY = (1.0, 1.0), dX = (1.0, 1.0)) #There is a slight error in the naming of X and Y variables
    nx, ny = size(u)
    #These are boundary conditions for all x's at the first position
    sum_dx = sum(dY)
    sum_dy = sum(dX)
    sum_d = sum_dx + sum_dy
    @inbounds for y in 2:ny-1
        x = 1
        du[x, y] += D * (
            sum_dx * u[x+1, y] + 
            dX[1] * u[x, y+1] + 
            dX[2] * u[x, y-1] - 
            sum_d * u[x, y])
    end

    #These are boundary conditions for all y's at the first position
    @inbounds for x in 2:nx-1
        y = 1
        du[x, y] += D * (
            dY[1] * u[x+1, y] + 
            dY[2] * u[x-1, y] + 
            sum_dy * u[x, y+1] - 
            sum_d * u[x, y])
    end

    #These are boundary conditions for x's at the end
    @inbounds for y in 2:ny-1
        x = nx
        du[x, y] += D * (
            sum_dx * u[x-1, y] + 
            dX[1] * u[x, y+1] + 
            dX[2] * u[x, y-1] - 
            sum_d * u[x, y])
    end

    #These are boundary conditions for all y's at end position
    @inbounds for x in 2:ny-1
        y = ny
        du[x, y] += D * (
            dY[1] * u[x+1, y] + 
            dY[2] * u[x-1, y] + 
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
            dY[1] * u[x+1, y] + 
            dY[2] * u[x-1, y] + 
            dX[1] * u[x, y-1] + 
            dX[2] * u[x, y+1] - 
            sum_d * u[x, y])
    end
end