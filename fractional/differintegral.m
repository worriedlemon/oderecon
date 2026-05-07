function Dx = differintegral(x, t, alpha, type)
    negative_alpha = alpha < 0;

    if (~exist('type', 'var'))
        type = 'gl';
    end

    Dx = zeros(size(x));
    Dx(:, negative_alpha) = int_rl(-alpha(negative_alpha), x(:, negative_alpha), t);
    Dx(:, ~negative_alpha) = fractdiff(x(:, ~negative_alpha), t, alpha(~negative_alpha), type);
end