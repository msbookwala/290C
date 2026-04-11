function ddx = ddx_bwd(f, dx)
    f = [f(1,:);f];
    ddx = (f(2:end,:)-f(1:end-1,:))/dx;
end