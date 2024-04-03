function orthogonality_test(F, deg, vc, a, b, eps)
    if ~exist('eps', 'var')
        eps = 1e-6;
    end
    
    N = size(F, 1);
    sigma2 = deglexord(deg * 2, vc);
    interv = repmat([a; b], 1, vc);
    
    R = zeros(N);
    for i = 1:N
        for j = 1:N
            R(i, j) = scalarpoly(F(i, :), F(j, :), sigma2, interv);
            if abs(R(i, j)) < eps
                R(i, j) = 0;
            end
        end
    end
    disp("Orthogonality test:");
    disp(R);
    disp("\nIf matrix is identity, then the transormation is completed successfully\n")
end