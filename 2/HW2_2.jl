using MAT
using Plots
using ColorSchemes
data = matread("2/gotritons_2.mat")
T = data["T"]'
x = data["xx"][:,1]
y = data["yy"][1,:]
# plot initial condition
# plot1 = heatmap(x, y, T, title="Temperature at t=0", xlabel="x", ylabel="y", colorbar_title="T", aspect_ratio=1);
# plot1

function ftbs_advection(T_init, dt, dx, dy, nt, alpha)
    # forward time, backward space
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

        T[n+1,2:end-1,2:end-1] = alpha[1]*dt/dx * T[n,1:end-2,2:end-1] + 
                                alpha[2]*dt/dy * T[n,2:end-1,1:end-2] +
                                (1 - alpha[1]*dt/dx - alpha[2]*dt/dy)* T[n,2:end-1,2:end-1]
    end
    return T[:,2:end-1,2:end-1]
end

function ftcs_advection(T_init, dt, dx, dy, nt, alpha)
    # forward time, backward space
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

        T[n+1,2:end-1,2:end-1] = -alpha[1]*dt/dx * (T[n,3:end,2:end-1] - T[n,1:end-2,2:end-1])/2 - 
                                alpha[2]*dt/dy * (T[n,2:end-1,3:end] - T[n,2:end-1,1:end-2])/2 +
                                T[n,2:end-1,2:end-1]
    end
    return T[:,2:end-1,2:end-1]
end

alpha = [1, 1]
final_t = 2
dx = x[2] - x[1]
dy = y[2] - y[1]
SF = 1.1
dt  = 1/(alpha[1]/dx + alpha[2]/dy)/SF
nt = ceil(Int, final_t/dt)
T_ftbs_ad = ftbs_advection(T, dt, dx, dy, nt, alpha)
# plot final condition
plot2 = heatmap(x, y, T_ftbs_ad[end,:,:], title="Temperature at t=$final_t",
                 xlabel="x", ylabel="y", colorbar_title="T",
                 seriescolor  = :viridis, aspect_ratio=1);
savefig(plot2, "2/advection_ftbs.png")



final_t = 0.25
nt = ceil(Int, final_t/dt)
T_ftcs_ad = ftcs_advection(T, dt, dx, dy, nt, alpha)
# plot final condition
plot3 = heatmap(x, y, T_ftcs_ad[end,:,:], title="Temperature at t=0.25",
                 xlabel="x", ylabel="y", colorbar_title="T",
                 seriescolor  = :viridis, aspect_ratio=1, zlims=(0,1));
savefig(plot3, "2/advection_ftcs.png")