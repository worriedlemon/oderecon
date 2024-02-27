function base = bernbase(t, sigma)
  [ic, vc] = size(sigma);
  vc = vc / 2;
  base = ones(size(t, 2), ic);
  for k = 1:ic
    for j = 1:vc
      n = sigma(k, j);
      i = sigma(k, j + vc);
      tx = t(j,:);
      base(:,k) = base(:,k).*(nchoosek(n, i).*tx.^i.*(1-tx).^(n-i))';
    end
  end
end