function p = EvalPolyOrth(h, X, T, F, sigma)
    [mc, ~] = size(F);
    h1 = zeros(mc, size(h, 2));
    for i = 1:mc
        flag = prod(sigma(i, :) == T, 2);
        if any(flag)
           [~, idn] = max(flag);
           h1(i, :) = h(idn, :);
        end
    end

    p = EvalPoly(F' * h1, X, sigma);
end