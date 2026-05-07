function Irv = int_rl(alpha, f, t, varargin)
    % INT_RL - Riemann-Liouville integral for alpha > 0
    %
    % - Irv = INT_RL(alpha, f, t)
    % - Irv = INT_RL(alpha, f, t, C)

    assert(all(alpha > 0), 'Riemann-Liouville integral is defined for alpha > 0');
    t = t(:);
    f = f(:);

    [n, vc] = size(f);
    Irv = zeros(n, vc);

    if (isscalar(alpha))
        alpha = alpha * ones(1, vc);
    end

    C = zeros(1, vc);
    if (size(varargin, 2) > 0)
        C = varargin{1, 1};
    end

    smallalpha = alpha < 1;
    bigalpha = ~smallalpha;

    for i = 2:n
        if (any(bigalpha))
            I = f(1:i, bigalpha) .* (t(i) - t(1:i)).^(alpha(bigalpha) - 1);
            Irv(i, bigalpha) = trapz(t(1:i), I);
        end

        if (any(smallalpha))
            I = ((t(i) - t(1:i-1)).^alpha(smallalpha) - (t(i) - t(2:i)).^alpha(smallalpha)) ./ alpha(smallalpha);
            Irv(i, smallalpha) = sum(I .* (f(1:i-1, smallalpha) + f(2:i, smallalpha)) / 2);
        end
    end

    Irv = Irv ./ gamma(alpha) + C;
end