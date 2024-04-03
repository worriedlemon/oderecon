function sc = scalarpoly(f, g, sigma2, interv)
    assert(size(interv) == [2, size(sigma2, 2)], 'Interval should be of size 2 x N, where N is variable count');
    
    a = interv(1, :); b = interv(2, :);
    cnt1 = length(f);
    cnt2 = length(g);
    
    sc = 0;
    for fmon = 1:cnt1
        if f(fmon) ~= 0
            pr = 1;
            shft = sigma2(fmon, :);
            for gmon = 1:cnt2
                if g(gmon) ~= 0
                    [~, ind] = max(prod(sigma2(gmon, :) + shft == sigma2, 2));
                    degs = sigma2(ind, :) + 1;
                    sc = sc + f(fmon) * g(gmon) * prod((b.^degs - a.^degs) ./ degs);
                end
            end
        end
    end
end