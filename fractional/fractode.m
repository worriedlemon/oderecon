function [t, y] = fractode(fun, t, y0, alpha)
    N = length(t);
    vc = length(y0);
    y = zeros(N, vc);
    y(1, :) = y0;
    
    w = zeros(1, N);
    w(1) = 1;
    for j = 2:N
        w(j) = -w(j-1) * (alpha - j + 2) / (j - 1);
    end

    for i = 2:N
        dy = zeros(1, vc);
        for j = 1:i-1
            s = w(j) * y(i - j, :);
            if (norm(s) < 1e-2)
                break;
            end
            dy = dy + s;
        end
        val = fun(t, y(i - 1, :));
        y(i, :) = (t(i) - t(i-1))^alpha * val' - dy;
    end
end
