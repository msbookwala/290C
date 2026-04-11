using MAT
using Plots
data = matread("1/cylinder_Re100.mat")
x = data["x"]
y = data["y"]
u = data["u"]
v = data["v"]
dt = data["dt"]
x_ = x[:,1]
y_ = y[1,:]

# part 1 #####
for i in 1:100
    plot1 = heatmap(x_, y_, u[1,:,:], title="Velocity in x-direction", xlabel="x", ylabel="y", colorbar_title="u", aspect_ratio=1)
    plot2 = heatmap(x_, y_, v[1,:,:], title="Velocity in y-direction", xlabel="x", ylabel="y", colorbar_title="v", aspect_ratio=1)
    plot!(plot1, plot2, layout=(2, 1), size=(1000, 1000))
end

