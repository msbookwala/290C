using MAT
using Plots
using ColorSchemes
data = matread("2/gotritons.mat")
T = data["T"]'
x = data["xx"][:,1]
y = data["yy"][1,:]
# plot initial condition
# plot1 = heatmap(x, y, T, title="Temperature at t=0", xlabel="x", ylabel="y", colorbar_title="T", aspect_ratio=1);
# plot1

function ftcs_diffusion(T_init, dt, dx, dy, nt, alpha)
    Tx_size = size(T_init) .+2
    T = zeros(nt+1, Tx_size...)
    T[1,2:end-1,2:end-1] = T_init
    for n in 1:nt
        #set periodic boundary conditions
        T[n,2:end-1,1] = T[n,2:end-1,end-1] # left boundary
        T[n,2:end-1,end] = T[n,2:end-1,2] # right boundary
        T[n,1,2:end-1] = T[n,end-1,2:end-1] # bottom boundary
        T[n,end,2:end-1] = T[n,2,2:end-1] # top boundary

        # 4 corner points
        T[n,1,1] = T[n,end-1,end-1] # bottom left corner
        T[n,1,end] = T[n,end-1,2] # bottom right corner
        T[n,end,1] = T[n,2,end-1] # top left corner
        T[n,end,end] = T[n,2,2] # top right corner
        # update interior points

        T[n+1,2:end-1,2:end-1] = alpha*dt*( T[n,1:end-2,2:end-1]/dx^2 + T[n,3:end,2:end-1]/dx^2 +
                                            T[n,2:end-1,1:end-2]/dy^2 + T[n,2:end-1,3:end]/dy^2 ) +
                                            T[n,2:end-1,2:end-1]*(1-2*alpha*dt/dx^2 - 2*alpha*dt/dy^2) 
    end

    return T[:,2:end-1,2:end-1]
end

alpha = 2
final_t = 0.001
dx = x[2] - x[1]
dy = y[2] - y[1]
SF = 1.001
dt  = 0.5/alpha/(1/dx^2 + 1/dy^2)/SF
nt = ceil(Int, final_t/dt)
T = ftcs_diffusion(T, dt, dx, dy, nt, alpha)
# plot final condition
plot2 = heatmap(x, y, T[end,:,:], title="Temperature at t=0.001",
                 xlabel="x", ylabel="y", colorbar_title="T",
                  aspect_ratio=1);
plot2