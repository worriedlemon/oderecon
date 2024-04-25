function coeff_restored = from01(coeff, domain, sigma, F)
    [mc, eqs] = size(coeff);
    [~, vc] = size(domain);
    
    coeff = coeff ./ prod(([-1, 1] * domain) .^ sigma, 2);
    coeff_restored = zeros(mc, eqs);
    
    if ~exist('F', 'var')
        F = eye(mc);
    end
    
    for eq = 1:eqs
        H1 = coeff(:, eq);
        for mon = 1:mc
            T1 = sigma(mon, :);
            F1 = F(mon, :);
            H2 = H1(mon);
            Ttemp = zeros(size(H1));
            for mon2 = 1:mc
                F2 = F1(mon2);
                for v = 1:vc
                    T2 = T1(v);
                    c = bincoeff(T2, 0:T2) .* ((-domain(1, vc)) .^ (0:T2));
                    Ttemp(1:length(c)) = F2 * H2 * fliplr(c);
                end
            end
            coeff_restored(:, eq) = coeff_restored(:, eq) + Ttemp;
        end
    end
end