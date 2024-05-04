function [F, nrms] = orthpoly(deg, vc, a, b, varargin)
    % -- F = ORTHPOLY(deg, vc, a, b)
    % -- F = ORTHPOLY(deg, vc, a, b, eps)
    % -- F = ORTHPOLY(deg, vc, a, b, eps, nrm)
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
    
    
    eps = 1e-12;
    nrm = 1;
    
    narg = nargin - 4;
    if narg > 0
        eps = varargin{1, 1};
    end
    if narg > 1
        nrm = varargin{1,2};
    end
    
    sigma = deglexord(deg, vc);    
    [N, ~] = size(sigma);
    F = eye(N);
    nrms = ones(N, 1);
    interv = repmat([a; b], 1, vc);
    
    sigma = [sigma; deglexord(deg + 1, deg * 2, vc)];
    
    for i = 1:N
        g = F(i, :);
        for k = 1:i-1
            temp = scalarpoly(g, F(k, :), sigma, interv) .* F(k, :);
            if ~nrm
                temp = temp / scalarpoly(F(k, :), F(k,:), sigma, interv);
            end
            F(i, :) = F(i, :) - temp;
        end
        
        c = sqrt(scalarpoly(F(i, :), F(i, :), sigma, interv));
        if nrm
            F(i, :) = F(i, :) / c;
        else
            nrms(i) = c;
        end
    end
    
    F = F .* (abs(F) > eps);
end
