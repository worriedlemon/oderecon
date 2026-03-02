function Ht = orthpoly_int_sindy(t,x,sigma,lambda_orth)
    % SINDy algorithm with Fourier series approach and integral criterion (OrthPoly-Int-SINDY)
    % Usage:
    %    -- Ht = ORTHPOLY_INT_SINDY(t, x, sigma, lambda)
    % Returns:
    %    -- Ht - system coefficients
    % Parameters:
    %    -- t - time
    %    -- x - nonlinear system states
    %    -- sigma - degree-lexicographic order ideal
    %    -- lambda - sparsification parameter (tolerance)
    
    %x signal
    rx = x; 
    % monomial count
    mc = size(sigma, 1);
    % eq count = variable count
    eqc = size(rx, 2); 
    vc = eqc; 
    
    
    F0 = orthpoly_t_integration(sigma, t, rx); % orthogonal polynomial matrix of integrals
    
    H0 = zeros(mc, eqc);
    E0 = EvalPoly(F0', rx, sigma);
    for i = 1:eqc
        for j = 1:mc
            H0(j, i) = h_ji_byint(E0(:, j),t,rx(:, i));
        end
    end
    
    Ht3 = F0' * H0; %ordinary H form orthogonal
    
    for jj = 1:5
        %calculate K1
        K1 = zeros(eqc,1);
        for i = 1:eqc
            K1(i) = trapz(t,rx(:,i))/sum(F0' * H0(:,i))/(t(end)-t(1));
        end
        H0 = zeros(mc, eqc);

        for i = 1:eqc
            F0 = orthpoly_t_integration(sigma, t, rx, K1(i)); % orthogonal polynomial matrix of integrals
            E0 = EvalPoly(F0', rx, sigma);
            for j = 1:mc
                H0(j, i) = h_ji_byint(E0(:, j),t,rx(:, i));
            end
            Ht3(:,i) = F0' * H0(:,i); %ordinary H form orthogonal
        end
    end


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
    
            % obtain new orthogonal matrix
            F_temp = zeros(size(F0));
    
            %F_temp(biginds, biginds) = orthpoly_t(sigma_temp, t, rx);
            F_temp(biginds, biginds) = orthpoly_t_integration(sigma_temp, t, rx); % orthogonal polynomial matrix of integrals
    
            % get new E
            E0 = EvalPoly(F_temp', rx, sigma);
            for j = 1:mc
                %Ht3(j, ind) = intdiff4(rx(:, ind),E(:, j));
                Ht3(j, ind) = h_ji_byint(E0(:, j),t,rx(:, ind));
            end
            Ht3(:, ind) = F_temp' * Ht3(:, ind); %get regular monomials
        end
        k = k + 1;
    end
    Ht = Ht3;
end