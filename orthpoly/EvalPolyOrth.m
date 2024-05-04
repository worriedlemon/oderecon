function p = EvalPolyOrth(h, X, T, F, sigma)
    [N, M] = size(X);
    [L, Q] = size(h);

    p = zeros(N,Q);
    idx = zeros(1, L);
    for i = 1:L
        [~, idn] = max(prod(sigma == T(i, :), 2));
        idx(i) = idn;
    end
    for i = 1:L
        p = p + h(i, :) .* EvalPoly(F(idx(i), :)', X, sigma);
    end
end