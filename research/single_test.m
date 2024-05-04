warning off;

sys = @Rossler;
coefs = [...
    0, 0, 0.2;...
    0, 1, 0;...
    -1, 0.2, 0;...
    -1, 0, -5.7;...
    0, 0, 0;...
    0, 0, 0;...
    0, 0, 1;...
    0, 0, 0;...
    0, 0, 0;...
    0, 0, 0];

Tmax = 100;
h = 0.01;
[~, x] = ode45(sys, 0:h:Tmax, [4 -2 0]);

vc = size(x, 2);
deg = 2;
sigma = deglexord(deg, vc);

N = int32(12);
N1 = int32(size(sigma, 1));

eps = 0.1;
tol = 1e-5;

Kiter = 1000;
Koff = 3;
rng_i(Koff);

sx = 0;
sorth = 0;

polynomial_options = {'x', 'orth'};
F = orthpoly(deg, deglexord(deg * 2, vc), 0, 1, 1e-9, 0);

% singlify
F = single(F);
x = single(x);
%x = double(x); % keeping signal in double
y = transpose(sys(0, x'));

for K = 1:Kiter
    wh = waitbar(K/Kiter);

    ri = randi(size(y, 1), 1, N);
    rx = x(ri,:); ry = y(ri,:);
    
    for po = 1:2
        opt = {polynomial_options{1, po}};
        if strcmp(opt{1,1}, 'orth')
            opt = {opt{1,1}, F, sigma};
        end
        
        H = zeros(N1, vc);
        
        E = EvalPoly(eye(N1), rx, sigma, opt);
        h0 = (E'*E)\E'*ry;
        for i = 1:vc
            [H(:, i), ~] = delMinorTerms(rx, ry(:,i), sigma, tol, opt, h0(:,i), 0);
        end
        
        switch opt{1,1}
            case 'x'
                sx = sx + all(all(abs(H - coefs) <= eps));
            case 'orth'
                sorth = sorth + all(all(abs(F' * H - coefs) <= eps));
        end
    end
end

close(wh);
disp("Success rate:")
fprintf("x -> %0.3f%%\n", sx / Kiter * 100)
fprintf("orth -> %0.3f%%\n", sorth / Kiter * 100)