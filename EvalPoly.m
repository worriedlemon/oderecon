function p = EvalPoly(h, X, T, varargin)
    %  p = EVALPOLY(h, X, T)
    %  p = EVALPOLY(h, X, T, opt)
    %
    %  Evaluates polynomial in datapoints X of size N x M, M is
    %  dimension, N is number of data points made up of monomials
    %  p - values of polynomes in X, N x Q
    %  h - coefficients by terms t, where t[i] is a monomial of the
    %  corresponding order, h size is L x Q
    %  T - ordered monomials, structure depends on polynomial basis (see opt)
    %  opt - defines polynomial basis type ('x', 'orth', 'bernstein'); default is 'x' 
    %
    %  If opt is default ('x'):
    %    T - L x M ordered monomials (e.g., w.r.t. degree-lexicographic order) like this, for 2-dimensional data:
    %     1    x    y    x^2   xy  y^2  x^2y  xy^2  x^2y^2
    %    [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2]
    %
    %  If opt is 'orth':
    %    T - L x M ordered monomials:
    %     f_i (x,y),     i >= 0
    %    [0 0; 1 0; 0 1; 2 0; 1 1; 0 2; 2 1;  1 2;  2 2]
    %    Also needed relation matrix F of size L x K, where K is number of initial basis polynomials;
    %    Values packed using cell like {opt, F}
    %
    %  If opt is 'bernstein':
    %    T - L x 2*M ordered monomials like this, for 2-dimensional data (second degree):
    %     (1-x)^2    2(1-x)x      x^2     (1-x)(1-y)    (1-x)y       x(1-y)       xy       (1-y)^2     2(1-y)y      y^2
    %    [2 0 0 0;   2 0 1 0;   2 0 2 0;   1 1 0 0;     1 1 0 1;    1 1 1 0;    1 1 1 1;   0 2 0 0;    0 2 0 1;    0 2 0 2...]
    
    opt = 'x';
    narg = nargin - 3;
    if narg > 0
        opt = varargin{1, 1};
        try
            if strcmp(opt{1,1}, 'orth')
                if size(opt, 2) < 3
                    error('Relation matrix and order ideal are needed')
                end
                F = opt{1,2};
                sigma = opt{1,3};
            end
            opt = opt{1,1};
        catch ME
            if strcmp(opt, 'orth')
                error("If option is 'orth', then relation matrix is also needed. Values should be packed in cell like {'orth', F}.")
            end
        end
    end

    [N, M] = size(X);
    [L, Q] = size(h);

    p = zeros(N,Q); %evaluated polynomial

    switch opt
        case 'x'
            for i = 1:L
                p = p + h(i,:) .* prod(X.^repmat(T(i,:), N, 1), 2);
            end
        case 'orth'
            p = EvalPolyOrth(h, X, T, F, sigma);
        case 'bernstein'
            p = EvalPolyBern(h, X, T);
        otherwise
            error('No such polynomial implementation')
    end
end
