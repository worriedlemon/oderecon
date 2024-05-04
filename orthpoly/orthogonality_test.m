function orthogonality_test(F, deg, vc, a, b, varargin)
    % -- orthogonality_test(F, deg, vc, a, b)
    % -- orthogonality_test(F, deg, vc, a, b, output)
    % -- orthogonality_test(F, deg, vc, a, b, output, eps)
    % -- orthogonality_test(F, deg, vc, a, b, output, eps, nrm)
    %     Runs an orthogonality test for given relation matrix
    %     F, knowing degree deg, variable count vc and
    %     orthogonality interval [a; b].
    %
    %     F - relation matrix
    %     deg - polynomial degree
    %     vc - dimension (variable count)
    %     a, b - orthogonality interval [a; b]
    %     output - 'verbose' or 'brief'
    %     eps - parameter for skipping values < eps
    %     nrm - norms of polynomials (default are ones)
    
    narg = nargin - 5;
    
    output = 'verbose';
    eps = 1e-6;
    nrm = ones(size(F, 1), 1);
    
    if narg > 0
       output = varargin{1,1};
    end
    
    if narg > 1
        eps = varargin{1,2};
    end
    
    if narg > 2
        nrm = varargin{1,3};
    end
    
    N = size(F, 1);
    sigma2 = deglexord(deg * 2, vc);
    interv = repmat([a; b], 1, vc);
    
    R = zeros(N);
    for i = 1:N
        for j = 1:N
            R(i, j) = scalarpoly(F(i, :) / nrm(i), F(j, :) / nrm(j), sigma2, interv);
        end
    end
    disp("Orthogonality test:");
    switch output
        case 'verbose'
            disp(R .* (abs(R) >= eps));
            disp("\nIf matrix is diagonal then the transormation is completed successfully\n");
        otherwise
            if isdiag(R .* (abs(R) >= eps))
                disp('Success');
            else
                disp('Failure');
            end
    end
end