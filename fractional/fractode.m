function [t, Y] = fractode(fun, t, y0, alpha)
    % FRACTODE - Solving fractional differential equations (FDE) system
    % using Adams-Bashforth method
    %
    % [t, y] = FRACTODE(fun, t, y0, alpha)

    
    % [t, y] = FRACTODE(fun, t, y0, alpha, useRk4)
    
    t = t(:);
    t0 = t(1);
    T = t(end);
    N = length(t) - 1;
    
    h = (T - t0) / N;
    [icc, vc] = size(y0);
    
    Y = zeros(N+1, vc);
    Y(1, :) = y0(1, :);
    
    % weights
    w = (1:N)'.^alpha - (0:(N-1))'.^alpha;
    
    % coefficient
    if (isscalar(alpha))
        alpha = alpha * ones(1, vc);
    end

    if (~exist('useRk4', 'var'))
        useRk4 = false;
    end

    ga = h.^alpha ./ gamma(alpha+1);
    
    % array for right side values, as we need the whole history
    F = zeros(N, vc);

    % initial conditions calculation
    ic = ceil(alpha);
    icr = max(ic);
    d0 = zeros(icr, vc);
    d0(1:icc, 1:vc) = y0;
    for i = 2:icr
        d0(i, i > ic) = 0;
    end
    d0 = d0 ./ gamma(1:icr)';
    D = EvalPoly(d0, t, (0:(icr-1))');
    
    if (useRk4)
        %Y(1:4, :) = RK4(@(t_x,x) (t(1:3)-t_x), t(1:3), Y(1, :));
    else
        for n = 1:3
            if (n == 1)
                F(n, :) = fun(t0, y0')';
            elseif (n == 2)
                F(n, :) = 1/2 * (3 * fun(t(n), Y(n, :)')' - fun(t(n-1), Y(n-1, :)')');
            else
                F(n, :) = 1/12 * (23 * fun(t(n), Y(n, :)')' - 16 * fun(t(n-1), Y(n-1, :)')' + 5 * fun(t(n-2), Y(n-2, :)')');
            end

            Y(n+1, :) = D(1, :) + ga .* sum(w(n - (0:(n-1)), :) .* F(1:n, :), 1);
        end
    end

    for n = 4:N
        F(n, :) = 1/24 * (55 * fun(t(n), Y(n, :)')' - 59 * fun(t(n-1), Y(n-1, :)')' + 37 * fun(t(n-2), Y(n-2, :)')' - 9 * fun(t(n-3), Y(n-3, :)')');

        Y(n+1, :) = D(1, :) + ga .* sum(w(n - (0:(n-1)), :) .* F(1:n, :), 1);
    end
end