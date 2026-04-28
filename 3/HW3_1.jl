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

cv = 718
cp = 1005
R = cp - cv
function prim2cons(rho, u, v, T, cv)
    return permutedims(stack([rho, rho.*u, rho.*v, rho.*(cv.*T + (u.^2 + v.^2)/2)]),(3,1,2))
end
function cons2prim(U,R,cv)
    rho = U[1,:,:]
    u = U[2,:,:] ./ rho
    v = U[3,:,:] ./ rho
    Et = U[4,:,:]
    e = (2*Et ./ rho - (u.^2 + v.^2))
    T = e ./ cv
    p = rho.*R.*T
    return rho, u, v, T, p, e, Et
end
function sutherland(T; mu0 = 1.735e-5, S1 = 110.4)
    return mu0 * (T / 298.15).^(3/2) .* (298.15 + S1) ./ (T .+ S1)
end
cons = prim2cons(rho, u, v, T, cv)
mu = sutherland(T)
rho, u, v, T, p, e, Et = cons2prim(cons, R, cv)

rhoplot = heatmap(x, y, rho, title="Density", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
uplot = heatmap(x, y, rho.*u, title="x-flux", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
vplot = heatmap(x, y, rho.*v, title="y-flux", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
Tplot = heatmap(x, y, T, title="Temperature", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
pplot = heatmap(x, y, p, title="Pressure", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
eplot = heatmap(x, y, e, title="Internal Energy", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
Etplot = heatmap(x, y, Et, title="Total Energy", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize  = 4,xtickfontsize = 4, ytickfontsize = 4)
muplot = heatmap(x, y, mu, title="Viscosity", xlabel="x", ylabel="y", seriescolor  = :viridis, titlefontsize = 6, xguidefontsize  = 4, yguidefontsize   = 4,xtickfontsize = 4,
ytickfontsize=4)
finplot= plot(rhoplot,
uplot,
vplot,
Tplot,
pplot,
eplot,
Etplot,
muplot;
layout=(4,
2))
savefig(finplot,
"3/final_plot.pdf")