warning off;

Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.2 0 0 0 0 0 0 0; 0.2 0 0 -5.7 0 0 1 0 0 0]';
sys = @Rossler;
sysname = 'Rossler';
Href = Hrossler;

Tmax = 50;
hs = 1e-4*10.^(0:0.5:3);

vc = 3;
eqc = 3;
deg = 2;
sigma = deglexord(deg, vc);
mc = size(sigma, 1);

N = length(hs);

errt = zeros(1, N);
erro = errt;
for k = 1:N
    disp(k);
    
    t = single(0:hs(k):Tmax);
    [~, x] = ode45(sys, t, [4 -2 0]);
    %y = diff4(x, t);
    y = [diff(x) / h; (x(end, :) - x(end - 1, :)) / h];

    F = single(orthpoly_t(sigma, t, x));

    E = single(EvalPoly(F', x, sigma));
    Ho = zeros(mc, eqc);
    for i = 1:eqc
        for j = 1:mc
            Ho(j, i) = trapz(x(:, i), E(:, j));
        end
    end

    Ho = F' * Ho;

    E = EvalPoly(eye(mc), x, sigma);
    Ht = (E'*E)\E'*y;

    errt(k) = norm(Ht - Href);
    erro(k) = norm(Ho - Href);
end

figure(1)
loglog(hs, errt + 1e-14, 'r', hs, erro + 1e-14, 'b');
grid on;
title(['Reconstruction with reduced precision (', sysname, ')']);
legend('LSM', 'Orthogonal Polynomials');
xlabel('Simulation time step \it{h}'); ylabel('Coefficients error \zeta');