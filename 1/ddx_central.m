function ddx = ddx_central(f, dx)
    f = [f(1,:);f;f(end,:)];
    ddx = (f(3:end,:)-f(1:end-2,:))/2/dx;
end