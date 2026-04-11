function ddx = ddx_fwd(f, dx)
    f = [f;f(end,:)];
    ddx = (f(2:end,:)-f(1:end-1,:))/dx;
end