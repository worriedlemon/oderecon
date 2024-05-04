function [F, nrms] = orthpoly_t(deg, vc, t, x_t, nrm)
    % -- F = ORTHPOLY(deg, vc, t, x_y)
    % -- F = ORTHPOLY(deg, vc, t, x_y, nrm)
    % -- [F, nrms] = ORTHPOLY(____)
    %     Returns relations matrix F, where the row
    %     describes polynomial coefficients. Also can
    %     return norms of every polynomial nrms.
    %
    %     deg - degree of a orthogonal polynomials
    %     vc - dimension (variables count)
    %     a, b - orthogonality interval edge values [a; b]
    %     eps - parameter for skipping values < eps
    %     nrm - logical value, indicating whether
    %       polynomials should be normalized (if so,
    %       dot product will be Kroneckers symbol)
    
    
    if ~exist('nrm', 'var')
        nrm = 1;
    end
    
    sigma = deglexord(deg, vc);
    N = size(sigma, 1);
    F = eye(N);
    nrms = ones(N, 1);
    
    for i = 1:N
        for k = 1:i-1
            n = EvalPoly(F(k, :)', x_t, sigma);
            temp = trapz(t, EvalPoly(F(i, :)', x_t, sigma) .* n) .* F(k, :);
            if ~nrm
                temp = temp / trapz(t, n .^ 2);
            end
            F(i, :) = F(i, :) - temp;
        end
        
        c = sqrt(trapz(t, EvalPoly(F(i, :)', x_t, sigma) .^ 2));
        if nrm
            F(i, :) = F(i, :) / c;
        else
            nrms(i) = c;
        end
    end
end
