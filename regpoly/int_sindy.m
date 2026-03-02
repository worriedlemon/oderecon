function Ht = int_sindy(t,x,sigma,lambda)
    % SINDy algorithm with integral criterion (Int-SINDY)
    % Usage:
    %    -- Ht = INT_SINDY(t, x, sigma, lambda)
    % Returns:
    %    -- Ht - system coefficients
    % Parameters:
    %    -- t - time
    %    -- x - nonlinear system states
    %    -- sigma - degree-lexicographic order ideal
    %    -- lambda - sparsification parameter (tolerance)
    
    [N, eqc] = size(x); % eq count = variable count
    %x signal
    rx = x; 
    
    % monomial count
    mc = size(sigma, 1);
    
    I = eye(mc);
    gam = zeros(mc,1);

    E_int = zeros(N,mc);
    %calculate mean of monomials
    for j = 1:mc
        E_int(:,j) = integrate_trapz(EvalPoly(I(:,j),rx,sigma),0,t);
        gam(j) = trapz(t, E_int(:,j)) / (t(end) - t(1));
    end
    
    Ht = zeros(mc, eqc);
    for i = 1:eqc
        x_avg = trapz(t,rx(:,i))/(t(end)-t(1));
        %V = rx(:,i) - rx(1,i);
        V = rx(:,i) - rx(1,i) - x_avg;
        Htemp = (E_int - gam')\V; %LSM
         % Original SINDy sparsification
        for iter = 1:10
            small = abs(Htemp) < lambda;
            Htemp(small) = 0;
            big = ~small;
            if sum(big) > 0
                Htemp(big) = (E_int(:,big) - (gam(big))') \ V;
            end
        end
        Ht(:,i) = Htemp;
    end
    
end