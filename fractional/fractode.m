function [t, Y] = fractode(fun, t, y0, alpha)
    % FRACTODE - Solving fractional differential equations (FDE) system
    % using Adams-Bashforth method
    
    t = t(:);
    t0 = t(1);
    T = t(end);
    N = length(t) - 1;
    
    h = (T - t0) / N;
    dim = length(y0);
    
    Y = zeros(N+1, dim);
    Y(1, :) = y0;
    
    % weights
    w = (1:N).^alpha - (0:(N-1)).^alpha;
    
    % coefficient
    ga = h^alpha / gamma(alpha+1);
    
    % array for right side values, as we need the whole history
    F = zeros(N, dim);
    
    for n = 1:N
        if (n == 1)
            % Euler
            F(1, :) = fun(t0, y0')';
        elseif (n == 2)
            % Adams-Bashfort 2 order
            F(n, :) = 3/2 * fun(t(n), Y(n, :)')' - 1/2 * fun(t(n-1), Y(n-1, :)')';
        else
            % Adams-Bashfort 3 order
            F(n, :) = 23/12 * fun(t(n), Y(n, :)')' - 16/12 * fun(t(n-1), Y(n-1, :)')' + 5/12 * fun(t(n-2), Y(n-2, :)')';
        end

        % Next point
        Y(n+1, :) = Y(1, :) + ga * w(n - (0:(n-1))) * F(1:n, :);
    end
end