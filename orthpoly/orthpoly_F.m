function [F, nn] = orthpoly_F(sigma, t, x_t, F, nstart, nn)
    % -- F = orthpoly_F(sigma, t, x_y, F, nstart)
    % -- [F, nn] = ORTHPOLY(____)
    %     Returns relations matrix F, where the row
    %     describes polynomial coefficients. This
    %     algorithm constructs orthogonal polynomials
    %     based on a time-dependent systems. Also can
    %     return norms of every polynomial nrms.
    %
    %     sigma - order ideal for polynomials
    %     F - matrix for orthogonalization
        
    N = size(sigma, 1);

    if ~exist('nn', 'var')
        nn = zeros(length(x_t),N);
        for i = 1:N
            nn(:,i) = EvalPoly(F(i, :)', x_t, sigma);
        end
    else
       for i = nstart:N
          nn(:,i) = EvalPoly(F(i, :)', x_t, sigma);
       end
    end
   
    for i = nstart:N %starting from a certain line in F
        for k = 1:i-1
            temp = trapz(t, nn(:,i) .* nn(:,k)) .* F(k, :);
            F(i, :) = F(i, :) - temp;
        end
        nn(:, i) = EvalPoly(F(i, :)', x_t, sigma);
        c = sqrt(trapz(t, nn(:, i) .^ 2));
        F(i, :) = F(i, :) / c;
        nn(:, i) = nn(:, i) / c;
    end
end
