function [F, nrms] = orthpoly_t_noint(sigma, x_t, nrm)
    
    if ~exist('nrm', 'var')
        nrm = 1;
    end
    
    N = size(sigma, 1);
    F = eye(N);
    nrms = ones(N, 1);
    
    for i = 1:N
        for k = 1:i-1
            n = EvalPoly(F(k, :)', x_t, sigma);
            E = EvalPoly(F(i, :)', x_t, sigma) .* n;
            temp = sum(E) .* F(k, :);
            if ~nrm
                temp = temp / sum(n .^ 2);
            end
            F(i, :) = F(i, :) - temp;
        end
        
        c = sqrt(sum(EvalPoly(F(i, :)', x_t, sigma) .^ 2));
        if nrm
            F(i, :) = F(i, :) / c;
        else
            nrms(i) = c;
        end
    end
end
