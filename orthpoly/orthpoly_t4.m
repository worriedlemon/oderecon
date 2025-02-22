function [F, nrms] = orthpoly_t4(sigma, t, x_t, nrm)
    % -- F = orthpoly_t4(sigma, t, x_y)
    % -- F = orthpoly_t4(sigma, t, x_y, nrm)
    % -- [F, nrms] = ORTHPOLY_T4(____)
    %     Returns relations matrix F, where the row
    %     describes polynomial coefficients. This
    %     algorithm constructs orthogonal polynomials
    %     based on a time-dependent systems. Uses 4 order
    %     integration method.
    %     Also can return norms of every polynomial nrms.
    %
    %     sigma - order ideal for polynomials
    %     a, b - orthogonality interval edge values [a; b]
    %     nrm - logical value, indicating whether
    %       polynomials should be normalized (if so,
    %       dot product will be Kroneckers symbol)
    
    
    if ~exist('nrm', 'var')
        nrm = 1;
    end
    
    N = size(sigma, 1);
    F = eye(N);
    nrms = ones(N, 1);
    
    for i = 1:N
        for k = 1:i-1
            n = EvalPoly(F(k, :)', x_t, sigma);
            temp = intdiff4(t, EvalPoly(F(i, :)', x_t, sigma) .* n) .* F(k, :);
            if ~nrm
                temp = temp / trapz(t, n .^ 2);
            end
            F(i, :) = F(i, :) - temp;
        end
        intgr = intdiff4(t, EvalPoly(F(i, :)', x_t, sigma) .^ 2);
        c = sqrt(intgr);
        if nrm
            F(i, :) = F(i, :) / c;
        else
            nrms(i) = c;
        end
    end
end
