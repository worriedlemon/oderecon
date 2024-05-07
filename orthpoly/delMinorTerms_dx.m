function [h, T] = delMinorTerms_dx(dx, x_t, ref_y, F, sigma, eta, deleteminor)

if ~exist('deleteminor', 'var')
    deleteminor = 1;
end

[N, ~] = size(x_t);
[L, M] = size(sigma);

h = zeros(L, 1);
for j = 1:L
    h(j) = trapz(dx, EvalPoly(F(j, :)', x_t, sigma));
end



T = sigma;
Ftmp = F;
F0 = F;
htmp = h;
Ttmp = T;

zeroingT = ones(L,1);
reindex = 1:L;
L0 = L;

while 1/N*norm(ref_y - EvalPolyOrth(htmp, x_t, Ttmp, F0, sigma)) <= eta && L > 1
    h = htmp;
    T = Ttmp;
    
    %find a minimal norm of the monomial in interpolation
    minval = inf;
    mink = 1;
    for k = 1:L
        t = zeros(L,1);
        t(k) = h(k); %extract k-th monomial
        tmp = norm(EvalPolyOrth(t, x_t, T, F0, sigma));
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
            Ftmp = Ftmp([1:mink - 1,mink + 1:end],:);
            reindex = reindex([1:mink - 1,mink + 1:end]);
        end
        if mink == L
            Ttmp = T(1:mink - 1,:);
            Ftmp = Ftmp(1:mink - 1,:);
            reindex = reindex(1:mink - 1);
        end
        if mink == 1
            Ttmp = T(mink + 1:end,:);
            Ftmp = Ftmp(mink + 1:end,:);
            reindex = reindex(mink + 1:end);
        end
        L = L - 1;
       
    end
    
    htmp = zeros(L, 1);
    for j = 1:L
        htmp(j) = trapz(dx, EvalPoly(Ftmp(j, :)', x_t, sigma));
    end
end

if 1/N*norm(ref_y - EvalPolyOrth(htmp, x_t, Ttmp, F0, sigma)) <= eta
    h = htmp;
    T = Ttmp;
elseif L < L0
    zeroingT(lmink) = 1;
end

if ~deleteminor
    htmp = zeros(L0,1);
    ctr = 1;
    for i = 1:L0
        if( zeroingT(i) ~= 0) && (ctr <= length(h))
            htmp(i) = h(ctr);
            ctr = ctr + 1;
        end
    end
    h = htmp;
    T = sigma;
end
end