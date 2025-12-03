function dX = VDPL(t,X)
    dX = X;
    mu = 1;
    
    x = X(1, :);
    y = X(2, :);
    
    dX(1, :) = y;
    dX(2, :) = mu * (1 - x.^2) .* y - x;
end