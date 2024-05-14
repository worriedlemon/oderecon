warning off;

Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.2 0 0 0 0 0 0 0; 0.2 0 0 -5.7 0 0 1 0 0 0]';
sys = @Lorenz;
sysname = 'Lorenz';
Href = Hlorenz;

Tmax = 50;
h = 1e-3;

vc = 3;
eqc = 3;
deg = 2;
sigma = deglexord(deg, vc);
mc = size(sigma, 1);

Kiter = 100;

errt = zeros(1, Kiter);
erro = errt;
for K = 1:Kiter
    disp(K);

    [t, x] = ode45(sys, 0:h:Tmax, randn(1, 3));
    x = single(x);
    y = diff4(x, t);
    
    F = orthpoly_t(sigma, t, x);

    E = EvalPoly(F', x, sigma);
    Ho = zeros(mc, eqc);
    for i = 1:eqc
        for j = 1:mc
            Ho(j, i) = trapz(x(:, i), E(:, j));
        end
    end

    Ho = F' * Ho;

    E = EvalPoly(eye(mc), x, sigma);
    Ht = (E'*E)\E'*y;

    errt(K) = norm(Ht - Href);
    erro(K) = norm(Ho - Href);
end

figure
subplot(121);
histogram(errt, 10, 'FaceColor', 'b');
grid on;
title(['Errors of reconstruction (', sysname, ')', newline, 'Method: LSM']);
ylabel('Norm error in coefficients');
subplot(122)
histogram(erro, 10, 'FaceColor', 'r');
grid on;
title(['Errors of reconstruction (', sysname, ')', newline, 'Method: Orthogonal Polynomials']);
ylabel('Norm error in coefficients');