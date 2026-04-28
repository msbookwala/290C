using MAT
using Plots
using ColorSchemes
data = matread("3/supersonicJetLES_xyPlane.mat")

u = data["u"]'
v = data["v"]'
x = data["xx"][:,1]
y = data["yy"][1,:]
xx = data["xx"]
yy = data["yy"]
rho = data["rho"]'
T = data["T"]'

cv = 1.0
function prim2cons(rho, u, v, T, cv)
    return permutedims(stack([rho, rho.*u, rho.*v, rho.*(cv.*T + (u.^2 + v.^2)/2)]),(3,1,2))
end

out = prim2cons(rho, u, v, T, cv)