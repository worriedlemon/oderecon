function F = orthpoly_t(sigma, t, x_t)
    % -- F = orthpoly_t(sigma, t, x_y)
    %     Returns relations matrix F, where the row
    %     describes polynomial coefficients. This
    %     algorithm constructs orthogonal polynomials
    %     based on a time-dependent systems. Also can
    %     return norms of every polynomial nrms.
    %
    %     sigma - order ideal for polynomials
    %     F - matrix for orthogonalization
        
    [F, ~] = orthpoly_F(sigma, t, x_t, eye(size(sigma, 1)), 1);
end
