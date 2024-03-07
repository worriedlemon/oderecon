function sigma = berndeg(n, vc)
  ns = deglexord(n, n, vc);
  
  sigma = [];
  for k = 1:size(ns, 1)
    nsc = ns(k,:);
    is = zeros(1, vc);
    loop = true;
    while loop
      sigma = [sigma; nsc is];
      try
        is = add_one(is, nsc);
      catch ME
        loop = ~loop;
      end
    end
  end
end

function x = add_one(x, s)
  last = length(x);
  x(last) = x(last) + 1;
  while ~prod(x <= s)
    [~, cur] = max(fliplr(x > s));
    cur = last - cur + 1;
    if (cur > 1)
      x(cur) = 0;
      x(cur - 1) = x(cur - 1) + 1;
    else
      error('Not possible to add another one')
    endif
  end
end