rng_i default;
warning off;

Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.2 0 0 0 0 0 0 0; 0.2 0 0 -5.7 0 0 1 0 0 0]';
system = @Lorenz;
sysname = 'Lorenz';
Href = Hlorenz;

Tmax = 100; % Time end
h = 1e-3; % Step
start_point = [0.1 0 0.1]; % Initial point

[t, x] = ode45(system, 0:h:Tmax, start_point);
y = transpose(system(0, x'));

[N, vc] = size(x);

eqc = size(y, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

noises = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.25, 0.5, 1, 2, 5, 10];
errt = []; erro = [];
for noise_amp = noises
    dmx = x - mean(x);
    rx = x + noise_amp * randn(N, vc) .* dmx;
    %ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
    %ry = diff2(rx) / h; % second order
    ry = diff4(rx, t); % fourth order
    
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    for i = 1:eqc
        E = EvalPoly(F', rx, sigma);
        for j = 1:mc
            %Ho(j, i) = trapz(t, E(:, j) .* ry(:, i));
            Ho(j, i) = trapz(rx(:, i), E(:, j));
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B)\B'*ry;
    %errt = [errt norm(y - EvalPoly(Ht, x, sigma))];
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];
    %erro = [erro norm(y - EvalPoly(Ht1, x, sigma))];
end

loglog(noises, errt, 'r', noises, erro, 'b');
grid on;
title(['Noise resistance (', sysname, ')']);
legend('LSM', 'Orthogonal polynomials');
xlabel('Noise magnitude');
ylabel('Norm error in coefficients');