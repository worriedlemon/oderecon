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
    w = (1:N)'.^alpha - (0:(N-1))'.^alpha;
    
    % coefficient
    ga = h.^alpha ./ gamma(alpha+1);
    
    % array for right side values, as we need the whole history
    F = zeros(N, dim);
    
    for n = 1:N
        if (n == 1)
            F(n, :) = fun(t0, y0')';
        elseif (n == 2)
            F(n, :) = 1/2 * (3 * fun(t(n), Y(n, :)')' - fun(t(n-1), Y(n-1, :)')');
        elseif (n == 3)
            F(n, :) = 1/12 * (23 * fun(t(n), Y(n, :)')' - 16 * fun(t(n-1), Y(n-1, :)')' + 5 * fun(t(n-2), Y(n-2, :)')');
        else
            F(n, :) = 1/24 * (55 * fun(t(n), Y(n, :)')' - 59 * fun(t(n-1), Y(n-1, :)')' + 37 * fun(t(n-2), Y(n-2, :)')' - 9 * fun(t(n-3), Y(n-3, :)')');
        end

        Y(n+1, :) = Y(1, :) + ga .* sum(w(n - (0:(n-1)), :) .* F(1:n, :), 1);
    end
end