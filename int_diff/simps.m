function intgr = simps(x, y)
    k = length(x);
    
    intgr = 0;
    for i = 1:2:k - 2
        Dx = (x(i) - x(i + 1)) * (x(i + 1) - x(i + 2)) * (x(i) - x(i + 2));
        a = (y(i) * (x(i + 1) - x(i + 2)) - y(i + 1) * (x(i) - x(i + 2)) + y(i + 2) * (x(i) - x(i + 1))) / Dx;
        b = -(y(i) * (x(i + 1)^2 - x(i + 2)^2) - y(i + 1) * (x(i)^2 - x(i + 2)^2) + y(i + 2) * (x(i)^2 - x(i + 1)^2)) / Dx;
        c = (y(i) * x(i + 1) * x(i + 2) * (x(i + 1) - x(i + 2)) - y(i + 1) * x(i) * x(i + 2) * (x(i) - x(i + 2)) + y(i + 2) * x(i) * x(i + 1) * (x(i) - x(i + 1))) / Dx;
        
        intgr = intgr + a / 3 * (x(i + 2) ^ 3 - x(i) ^ 3) + b / 2 * (x(i + 2) ^ 2 - x(i) ^ 2) + c * (x(i + 2) - x(i));
    end

    if i + 2 < k
        intgr = intgr + (x(k) - x(k - 1)) / 2 * (y(k - 1) + y(k));
    end
end