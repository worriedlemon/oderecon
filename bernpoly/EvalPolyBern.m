function p = EvalPolyBern(h, X, T)
    [N, M] = size(X);
    [L, Q] = size(h);

    p = zeros(N,Q); %evaluated polynomial
    for i = 1:L
        ns = T(i, 1:M);
        is = T(i, M + 1:2*M);
        mon = ones(N, Q);
        for k = 1:M
            mon = mon .* nchoosek(ns(k), is(k)).*X(:, k).^is(k).*(1 - X(:, k)).^(ns(k) - is(k));
        end
        p = p + h(i,:).*mon;
        % Possible faster formula
        % p = p + h(i, :) .* bincoeff(ns, is) .* (X .^ is) .* ((1 - X) .^ (ns - is));
    end
end