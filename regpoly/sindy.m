function H = sindy(x, y, sigma, lambda)
% Sparse Identification of Nonlinear Dynamics algorithm (SINDY)
% Usage:
%    -- H = SINDY(x, y, sigma, lambda)
% Returns:
%    -- H - system coefficients
% Parameters:
%    -- x - nonlinear system data
%    -- y - nonlinear system derivatives
%    -- sigma - degree-lexicographic order ideal
%    -- lambda - sparsification parameter (tolerance)

    delta = 0.01;
    [mc, vc] = size(sigma);
    B = EvalPoly(eye(mc), x, sigma);
    
    Theta = (B'*B + delta*eye(mc));
    D = B'*y;
    H = Theta \ D; % initial guess: Least-squares
    
    while (1)
        smallinds = (abs(H) < lambda); % minor indexes are going to be zeroed
        
        if (all(H(smallinds) == 0, "all")) % all values are zeros, stop
            break
        end

        H(smallinds) = 0;
        for ind = 1:vc
            biginds = ~smallinds(:, ind);
            H(biginds, ind) = Theta(:, biginds) \ D(:, ind); 
        end
    end
end