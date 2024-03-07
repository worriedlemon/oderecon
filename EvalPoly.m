function p = EvalPoly(h, X, T, opt)
    %  p = EVALPOLY(h, X, T) evaluates polynomial in datapoints X of size N x M, M is
    %  dimension, N is number of data points made up of monomials
    %  p - values of polynomes in X, N x Q
    %  h - coefficients by terms t, where t[i] is a monomial of the
    %  corresponding order, h size is L x Q
    %  T - ordered monomials, structure depends on polynomial basis (see opt)
    %  opt is optional parameter, defining polynomial basis type ('x' or 'bernstein'); default is 'x' 
    %
    %  If opt is default ('x'):
    %    T - L x M ordered monomials (e.g., w.r.t. degree-lexicographic order) like this, for 2-dimensional data:
    %     1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
    %    [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2]
    %  If opt is 'bernstein':
    %    T - L x 2*M ordered monomials like this, for 2-dimensional data (second degree):
    %     (1-x)^2    2(1-x)x      x^2     (1-x)(1-y)    (1-x)y       x(1-y)       xy       (1-y)^2     2(1-y)y      y^2
    %    [2 0 0 0;   2 0 1 0;   2 0 2 0;   1 1 0 0;     1 1 0 1;    1 1 1 0;    1 1 1 1;   0 2 0 0;    0 2 0 1;    0 2 0 2...]
    
    %get dimension
    if ~exist('opt', 'var')
        opt = 'x';
    end

    [N, M] = size(X);
    [L, Q] = size(h);

    p = zeros(N,Q); %evaluated polynomial

    switch opt
        case 'x'
            for i = 1:L
               p = p + h(i,:).*prod(X.^repmat(T(i,:),N,1),2);
            end
        case 'bernstein'
            for i = 1:L
               ns = T(i, 1:M);
               is = T(i, M + 1:2*M);
               mon = ones(N, Q);
               for k = 1:M
                 mon = mon .* nchoosek(ns(k), is(k)).*X(:, k).^is(k).*(1 - X(:, k)).^(ns(k) - is(k));
               end
               p = p + h(i,:).*mon;
            end
        otherwise
            error('No such polynomial implementation')
    end
end
