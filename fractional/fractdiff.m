function dx = fractdiff(x, t, alpha)
    % Grunwald-Letnikov differintegral

    assert(length(t) > 1, 't length should be > 1')

    [N, M] = size(x);
    
    % recurrent
    w = zeros(1, N);
    w(1) = 1;
    for j = 2:N
        w(j) = w(j-1) * (1 - (1 + alpha) / (j - 1));
    end
    
    dx = zeros(N, M);
    dx(1, :) = x(1, :) * (t(2) - t(1))^(-alpha);
    
    for i = 2:N
        h = (t(i) - t(i - 1));
        for j = 1:i
            s = w(j) * x(i - j + 1, :);
            if (norm(s) < h * h)
                break;
            end
            dx(i, :) = dx(i, :) + s;
        end
        
        dx(i, :) = dx(i, :) * h^(-alpha);
    end
end