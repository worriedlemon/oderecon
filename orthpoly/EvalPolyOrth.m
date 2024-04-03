function p = EvalPolyOrth(h, X, F, sigma)
    [N, M] = size(X);
    [L, Q] = size(h);

    p = zeros(N,Q);
    fs = EvalPoly(eye(L), X, sigma);
    for i = 1:L
        p = p + h(i, :) .* (F(i, :) .* fs);
    end
end