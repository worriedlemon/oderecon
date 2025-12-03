function Ht = orthpoly_sindy(t,x,sigma,lambda_orth)
    % Sparse Identification of Nonlinear Dynamics algorithm with Fourier series approach (OrthPoly-SINDY)
    % Usage:
    %    -- Ht = ORTHPOLY_SINDY(t, x, sigma, lambda)
    % Returns:
    %    -- Ht - system coefficients
    % Parameters:
    %    -- t - time
    %    -- x - nonlinear system states
    %    -- sigma - degree-lexicographic order ideal
    %    -- lambda - sparsification parameter (tolerance)

    rx = x; %x signal
    F = orthpoly_t(sigma, t, rx); % orthogonal polynomial matrix
    
    % monomial count
    mc = size(sigma, 1);
    % eq count = variable count
    eqc = size(x, 2); 
    vc = eqc; 
    
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc
        for j = 1:mc
            Ho(j, i) = intdiff4(rx(:, i), E(:, j));
        end
    end
    
    Ht3 = F' * Ho; %ordinary H form orthogonal
    k = 0;
    % avoiding infinite loop (probably, won't ever happen)
    while (k < 10)
        smallinds = (abs(Ht3) < lambda_orth);
        if (all(Ht3(smallinds) == 0, "all")) % all small values are already zeros, stop
            break
        end
    
        Ht3(smallinds) = 0;
        for ind = 1:vc
            biginds = ~smallinds(:, ind);
            sigma_temp = sigma(biginds, :);
    
            %obtain new orthogonal matrix
            F_temp = zeros(size(F));
            F_temp(biginds, biginds) = orthpoly_t(sigma_temp, t, rx);
    
            % get new E
            E = EvalPoly(F_temp', rx, sigma);
            for j = 1:mc
                Ht3(j, ind) = intdiff4(rx(:, ind),E(:, j));
            end
            Ht3(:, ind) = F_temp' * Ht3(:, ind); %get regular monomials
        end
        k = k + 1;
    end
    Ht = Ht3;
end