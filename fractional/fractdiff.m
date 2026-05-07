function dx = fractdiff(x, t, alpha, type)
    % FRACTDIFF - compute fractial derivative using differintegral
    %
    % - FRACTDIFF(x, t, alpha) - calculate fractional derivative d^alpha x/dt^alpha
    % - FRACTDIFF(x, t, alpha, type) - calculate fractional derivative
    %        d^alpha x/dt^alpha with specified derivative (default: 'gl')
    %
    % Arguments:
    %   x     - left side values (points) to take derivative from [N,M]
    %   t     - argument values to integrate along [N,1]
    %   alpha - differintegral order [1,M]
    %   type  - type of differintegral:
    %           * 'gl' for Grunwald-Letnikov
    %           * 'rl' for Riemann-Liouville
    %           * 'cp' for Caputo
    %           * 'cab' for Caputo/Adams-Bashforth

    assert(length(t) > 1, 't length should be > 1')
    [N, M] = size(x);

    if (~exist('type', 'var'))
        type = 'gl';
    end

    if (isscalar(alpha))
        alpha = alpha * ones(1, M);
    end

    nalpha = ceil(alpha);
    idx = 1:M;
    ints = abs(nalpha - alpha) < 1e-6;
    nintidx = idx(~ints);
    
    dx = zeros(N, M);
    
    % compute integer orders
    for j = idx(ints)
        dx(:, j) = fcompose(@(x)diff4(x, t), nalpha(j), x(:, j)); 
    end

    if (isempty(nintidx))
        return
    end
    
    % first point
    dx(1, nintidx) = gamma(alpha(nintidx) + 1) .* (t(2) - t(1)).^(-alpha(nintidx)) .* (x(2, nintidx) - x(1, nintidx));

    switch type
        case 'gl'
            % recurrent formula for weights
            w = ones(N, M);
            for j = 2:N
                w(j, :) = w(j-1, :) .* (1 - (1 + alpha) / (j - 1));
            end
            
            for i = 2:N
                idx = 1:i;
                dx(i, nintidx) = (t(i) - t(i - 1)).^(-alpha(nintidx)) .* sum(w(idx, nintidx) .* x(i + 1 - idx, nintidx), 1);
            end
        case 'rl'
            % using straight Riemann-Liouville formula
            dx(:, nintidx) = int_rl(nalpha - alpha, x(:, nintidx), t);
            for j = nintidx
                dx(:, j) = fcompose(@(x)diff4(x, t), nalpha(j), dx(:, j)); 
            end
        case 'cp'
            % using straight Caputo formula
            for j = nintidx
                dx(:, j) = fcompose(@(x)diff4(x, t), nalpha(j), x(:, j)); 
            end
            dx(:, nintidx) = int_rl(nalpha - alpha, dx(:, nintidx), t);
        case 'cab'
            % weights pre-calculation
            w = (1:N-1)'.^alpha - (0:N-2)'.^alpha;
    
            for i = 1 : N-1
                dx(i, nintidx) = gamma(alpha(nintidx) + 1) .* (t(i+1) - t(i)).^(-alpha(nintidx)) .* (x(i+1, nintidx) - x(1, nintidx));
                if (i > 1)
                    idx = 1:i-1;
                    dx(i, nintidx) = dx(i, nintidx) - sum(w(idx + 1, nintidx) .* dx(i - idx, nintidx), 1);
                end
            end
        otherwise
            if (ischar(type) || isstring(type))
                error(['Type ''', type, ''' is not supported']);
            else
                error('Type should be a string or character array');
            end
    end
end