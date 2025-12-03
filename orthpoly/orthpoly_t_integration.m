function [F, nrms] = orthpoly_t_integration(sigma, t, x_t, nrm)
    % -- F = orthpoly_t_integration(sigma, t, x_y)
    % -- F = orthpoly_t_integration(sigma, t, x_y, nrm)
    % -- [F, nrms] = ORTHPOLY_T_INTEGRATION(____)
    %     Returns relations matrix F, where the row
    %     describes polynomial coefficients. This
    %     algorithm constructs orthogonal polynomials
    %     based on a time-dependent systems. Also can
    %     return norms of every polynomial nrms.
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
            %n = EvalPoly(F(k, :)', x_t, sigma); 
            n = integrate_trapz(EvalPoly(F(k, :)', x_t, sigma),0,t); % integral of theta_k
            n1 = integrate_trapz(EvalPoly(F(i, :)', x_t, sigma),0,t);

            %n = detrend(n,0);
            n0 = trapz(t,n)/(t(end) - t(1));
            n = n - n0;

            %n1 = detrend(n1,0);
            n10 = trapz(t,n1)/(t(end) - t(1));
            n1 = n1 - n10;

            temp = trapz(t, n1 .* n) .* F(k, :);
            if ~nrm
                temp = temp / trapz(t, n .^ 2);
            end
            F(i, :) = F(i, :) - temp;
        end
        ni = integrate_trapz(EvalPoly(F(i, :)', x_t, sigma),0,t);
        %ni = detrend(ni,0);
        ni0 = trapz(t,ni)/(t(end) - t(1));
        ni = ni - ni0;

        intgr = trapz(t, ni .^ 2);
        c = sqrt(intgr);
        if nrm
            F(i, :) = F(i, :) / c;
        else
            nrms(i) = c;
        end

        %check for orthogonality
        % ni = integrate_trapz(EvalPoly(F(i, :)', x_t, sigma),0,t);
        % ni = detrend(ni,0);
        % for k = 1:i-1
        %     n = integrate_trapz(EvalPoly(F(k, :)', x_t, sigma),0,t); % integral of theta_k
        %     n = detrend(n,0);
        %     temp = trapz(t, ni .* n);
        %     if abs(temp) > 2e-11
        %         rrr = 0;
        %     end
        % end
    end
end
