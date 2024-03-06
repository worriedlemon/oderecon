function dx = oderecon(H,T,t,x,opt)
% dx = ODERECON(H,T,t,x)
% dx = ODERECON(H,T,t,x,opt)
%
%ODERECON uses coefficients h and set of basis polynmials T for
%reconstructing right-hand side of ODE
% x' = f(t,x)
% dx = ODERECON(H,T,t,x,opt) returns derivative vector x' in the point (t,x), where
%   H is a cell array of column vectors for polynomials H{1,i} by x'(i)
%   T is a cell array basis polyminials, with each row corresponding to
%   x'(i), e.g. T = [2 1] stands for  x(1)^2*x(2)
%   t is a scalar time value
%   x is a colunm vector
%   opt is an option for basis of evaluated polynomial (can be "x" or "bernstein"), default "x"
    if ~exist("opt")
        opt = "x";
    end
    dx = x;
    
    switch opt
        case "x"
            for i = 1:size(H, 2)
                dx(i) = EvalPoly(H{1,i}, x', T{1,i});
            end
        case "bernstein"
            for i = 1:size(H, 2)
                dx(i) = EvalBernstein(H{1,i}, x', T{1,i});
            end
        otherwise
            error("No such polynomial implementation");
    end
end
