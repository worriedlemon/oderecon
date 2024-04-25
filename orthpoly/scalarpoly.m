function sc = scalarpoly(f, g, sigma2, interv)
    % -- sc = scalarpoly(f, g, sigma2, interv)
    %     Returns dot product of two polynomials with
    %     weight function w(x) = 1 on a given interval.
    %
    %     f - polynomial of degree P1, dimension D1
    %     g - polynomial of degree P2, dimension D2
    %     sigma2 - order ideal of order `P1 + P2` and
    %       dimension `max(D1, D2)`
    %     interv - integration intervals in form of
    %       2 x max(D1, D2) size.
    
    assert(size(interv) == [2, size(sigma2, 2)], 'Interval should be of size 2 x N, where N is variable count');
    
    a = interv(1, :); b = interv(2, :);
    cnt1 = length(f);
    cnt2 = length(g);
    
    sc = 0;
    for fmon = 1:cnt1
        if f(fmon) ~= 0
            for gmon = 1:cnt2
                if g(gmon) ~= 0
                    [~, ind] = max(prod(sigma2(gmon, :) + sigma2(fmon, :) == sigma2, 2));
                    degs = sigma2(ind, :) + 1;
                    sc = sc + f(fmon) * g(gmon) * prod((b.^degs - a.^degs) ./ degs);
                end
            end
        end
    end
end