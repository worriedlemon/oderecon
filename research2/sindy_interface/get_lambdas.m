function lambdas = get_lambdas(methods, sys, sigma, Href)
    h = 1e-3;
    [t, x_ref] = ode78(sys, 0:h:100, load_pref(func2str(sys)));
    
    mtdn = size(methods, 2);
    lambdas = zeros(1, mtdn);
    test_count = 30;
    lv = logspace(-6, -1, test_count);
    vs = zeros(mtdn, test_count);
    figure(1);
    for k = 1:mtdn
        mtd = methods{1, k};
        disp(['Testing method ', method2str(mtd), '...']);
        optf = @(lambda)norm(Href - mtd(t, x_ref, sigma, lambda));
        for l = 1:test_count
            vs(k, l) = optf(lv(l));
        end
        [~, mi] = min(vs(k, :));
        lambdas(k) = lv(mi);
        disp(['Best lambda = ', num2str(lambdas(k))]);
    end
    semilogx(lv, vs);
    legend show;
end