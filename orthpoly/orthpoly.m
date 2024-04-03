function F = orthpoly(deg, vc, a, b)
    % Here the scalar product is equivalent to Kroneker delta function
    % (f_i, f_j) = 0 if i != j, otherwise 1
    
    sigma = deglexord(deg, vc);    
    [N, ~] = size(sigma);
    F = eye(N);
    interv = repmat([a; b], 1, vc);
    
    sigma = [sigma; deglexord(deg + 1, deg * 2, vc)];
    
    for i = 1:N
        g = F(i, :);
        for k = 1:i-1
            F(i, :) = F(i, :) - scalarpoly(g, F(k, :), sigma, interv) .* F(k, :);
        end
        F(i, :) = F(i, :) ./ sqrt(scalarpoly(F(i, :), F(i, :), sigma, interv));
    end
end
