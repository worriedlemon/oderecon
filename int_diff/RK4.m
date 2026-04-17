function [tspan,y] = RK4(fun,tspan,x0)

N = length(tspan);
dim = length(x0);
y = zeros(N,dim);
y(1,:) = x0;

b1 = 1/6;
b2 = 1/3;
b3 = 1/3;
b4 = 1/6;

x = x0;
for i = 2:N
    h = tspan(i) - tspan(i-1);
    
    t0 = tspan(i);
    k1 = fun(t0,x);
    x2 = x + 0.5*h*k1;
    k2 = fun(t0 + 0.5*h, x2);
    x3 = x + 0.5*h*k2;
    k3 = fun(t0 + 0.5*h, x3);
    x4 = x + h*k3;
    k4 = fun(t0 + h, x4);
    x = x + h*(b1*k1 + b2*k2 + b3*k3 + b4*k4);

    y(i,:) = x';
    if any(isnan(x))
        y = y(1:i-1,:);
        tspan = tspan(1:i-1);
        return;
    end
end

if (size(tspan, 1) == 1)
    tspan = tspan';
end

end