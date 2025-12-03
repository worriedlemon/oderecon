function dX = florenz(t,X,params)
    sig = params(1);
    bet = params(2);
    rho = params(3);
    dX = X;
    x = X(1); y = X(2); z = X(3);
    dX(1) = sig*(y - x);
    dX(2) = x*(rho - z) - y;
    dX(3) = x*y - bet*z;
end