rng_i default;
warning off;

Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.2 0 0 0 0 0 0 0; 0.2 0 0 -5.7 0 0 1 0 0 0]';
system = @Rossler;
sysname = 'Rossler';
Href = Hrossler;

Tmax = 50;
Tmaxs = 10:10:100; % Time end
hs = 10.^(-4:0.5:-1);
h = 1e-3; % Step
start_point = [0.1 0 0.1]; % Initial point

nrm = zeros(3, length(hs));
for i = 1:length(hs)
    disp(i)
    [t, x] = ode45(system, 0:hs(i):Tmax, start_point);
    y = transpose(system(0, x'));
    
    [N, vc] = size(x);
    
    eqc = size(y, 2);
    deg = 2;
    
    sigma = deglexord(deg, vc);
    mc = size(sigma, 1);
    
    B = EvalPoly(eye(mc), x, sigma);
    Ht = (B'*B)\B'*y;

    nrm(1, i) = 1 / (norm(Ht - Href) ^ 2 + 1);
    
    F = orthpoly_t(sigma, t, x);

    Ho_t = zeros(mc, eqc);
    Ho_x = Ho_t;
    for k = 1:eqc
        for j = 1:mc
            Ho_t(j, k) = trapz(t, EvalPoly(F(j, :)', x, sigma) .* y(:, k));
            Ho_x(j, k) = trapz(x(:, k), EvalPoly(F(j, :)', x, sigma));
        end
    end

    Ho_t = F' * Ho_t;
    Ho_x = F' * Ho_x;

    nrm(2, i) = 1 / (norm(Ho_t - Href) ^ 2 + 1);
    nrm(3, i) = 1 / (norm(Ho_x - Href) ^ 2 + 1);
end

figure(1);
%plot(Tmaxs, nrm(1, :), 'b--', Tmaxs, nrm(2, :), 'r--', Tmaxs, nrm(3, :), 'g');
semilogx(hs, nrm(1, :), 'b--', hs, nrm(2, :), 'r--', hs, nrm(3, :), 'g');
grid on;
ylim([0 1]);
title(['Reconstruction accuracy (', sysname, ')']);
%xlabel('Time end \it{T}_{max}, s');
xlabel('Time step \it\Delta{t}, s');
ylabel('Accuracy \it{\xi}');
legend('LSM', 'Orthogonal polynomials (time)', 'Orthogonal polynomials (var)')

figure(2);
%plot(Tmaxs, nrm(1, :), 'b--', Tmaxs, nrm(2, :), 'r--');
semilogx(hs, nrm(1, :), 'b--', hs, nrm(2, :), 'r--');
grid on;
title(['Reconstruction accuracy (', sysname, ')']);
%xlabel('Time end \it{T}_{max}, s');
xlabel('Time step \it\Delta{t}, s');
ylabel('Accuracy \it{\xi}');
legend('LSM', 'Orthogonal polynomials (time)');