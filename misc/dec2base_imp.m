function bch = dec2base_imp(X, K)
   if (K <= 36)
       bch = dec2base(X, K) - '0';
   else
       digitsN = floor(log(X) / log(K)) + 1;
       bch = zeros(1, digitsN);
       for i = 1:digitsN
           bch(digitsN - i + 1) = mod(X, K);
           X = floor(X / K);
       end
   end
end