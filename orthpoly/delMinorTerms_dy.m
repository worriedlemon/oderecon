function [h, T, F, h_reg] = delMinorTerms_dy(ts, dx, x_t, ref_y, F, sigma, eta, deleteminor)

if ~exist('deleteminor', 'var')
    deleteminor = 1;
end

[N, ~] = size(x_t);
[L, ~] = size(sigma);

h = zeros(L, 1);
for j = 1:L
    h(j) = trapz(dx, EvalPoly(F(j, :)', x_t, sigma));
end

T = sigma;
Ftmp = F;
htmp = h;
Ttmp = T;
h_reg = F'*h;

zeroingT = ones(L,1);
reindex = 1:L;
L0 = L;

while 1/N*norm(ref_y - EvalPoly(Ftmp'*htmp,x_t,Ttmp)) <= eta && L > 1
    F = Ftmp;
    h = htmp;
    h_reg = F'*h; %regular polynomials
    T = Ttmp;
    
    %find a minimal norm of the original (small) monomial in interpolation
    minval = inf;
    mink = 1;

    for k = 1:L
        t = zeros(L,1);
        t(k) = h_reg(k); %extract k-th monomial
        tmp = norm(EvalPoly(t,x_t,T));
        if tmp < minval
            minval = tmp;
            mink = k;
        end
    end
    
    lmink = reindex(mink);
    
    %exclude mink-th element from T
    if L > 1
        zeroingT(reindex(mink)) = 0;
        
        if mink > 1 && mink < L
            Ttmp = T([1:mink - 1,mink + 1:end],:);
            Ftmp = Ftmp([1:mink - 1,mink + 1:end],[1:mink - 1,mink + 1:end]);
            reindex = reindex([1:mink - 1,mink + 1:end]);
        end
        if mink == L
            Ttmp = T(1:mink - 1,:);
            Ftmp = Ftmp(1:mink - 1, 1:mink - 1);
            reindex = reindex(1:mink - 1);
        end
        if mink == 1
            Ttmp = T(mink + 1:end,:);
            Ftmp = Ftmp(mink + 1:end, mink + 1:end);
            reindex = reindex(mink + 1:end);
        end
        L = L - 1;
       
    end
    
    %Ftmp = orthpoly_t(Ttmp, ts, x_t); % Getting orthogonal polynomials matrix
    [Ftmp,~] = orthpoly_F(Ttmp, ts, x_t, Ftmp, mink); % Getting orthogonal polynomials matrix
    htmp = zeros(L,1);
    for j = 1:L
        htmp(j) = trapz(dx, EvalPoly(Ftmp(j, :)', x_t, Ttmp));
    end
end

if 1/N*norm(ref_y - EvalPoly(Ftmp'*htmp,x_t,Ttmp)) <= eta
    h = htmp;
    T = Ttmp;
    F = Ftmp;
    h_reg = F'*h;
elseif L < L0
    zeroingT(lmink) = 1;
end

if ~deleteminor
    htmp = zeros(L0,1);
    Ftmp = zeros(L0);
    ctr = 1;
    for i = 1:L0
        if( zeroingT(i) ~= 0) && (ctr <= length(h))
            htmp(i) = h(ctr);
            ctr2 = 1;
            for j = 1:L0
                if( zeroingT(j) ~= 0) && (ctr <= length(h))
                    Ftmp(i,j) = F(ctr,ctr2);
                    ctr2 = ctr2 + 1;
                end
            end
            ctr = ctr + 1;
        end
    end
    h = htmp;
    T = sigma;
    F = Ftmp;
    h_reg = F'*h;
end
end