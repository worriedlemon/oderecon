rng_i default;
warning off;

% System coefficients
Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.3 0 0 0 0 0 0 0; 0.3 0 0 -5.7 0 0 1 0 0 0]';

% Used system
system = @Rossler

Href = eval(['H', lower(func2str(system))]); % Used coefficients
Tmax = 100; % Time end
h = 1e-2; % Step
start_point = [4 -2 0]; % Initial point

[t, x] = ode45(system, 0:h:Tmax, start_point);
y = transpose(system(0, x'));

[N, vc] = size(x);

eqc = size(y, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

delta = 0.01; % regularization parameter

noises = logspace(-5, 1, 50);

% Homoscedasticity test
errt = []; erro = [];
for noise_amp = noises
    rx = x + noise_amp * randn(N, vc);
    ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
    %ry = diff2(rx) / h; % second order
    %ry = diff4(rx, t); % fourth order
    
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
    Ht = (B'*B + delta*eye(mc))\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];
end

figure(1);
loglog(noises, errt, 'r', noises, erro, 'b');
hold on; grid on;
%title(['Noise resistance (', func2str(system), ')']);
legend('LSM (Homoscedastic)', 'OrthPoly (Homoscedastic)');
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
%xlabel('Average noise magnitude $\sigma$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

% Heteroscedasticity test
errt = []; erro = [];
for noise_amp = noises
    dmx = x - mean(x);
    rx = x + noise_amp * randn(N, vc) .* dmx;
    ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h];
    
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    for i = 1:eqc
        E = EvalPoly(F', rx, sigma);
        for j = 1:mc
            Ho(j, i) = trapz(rx(:, i), E(:, j));
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B)\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];
end

figure(1);
hold on;
loglog(noises, errt, 'm--', 'DisplayName', 'LSM (Heteroscedastic)')
loglog(noises, erro, 'g--', 'DisplayName', 'OrthPoly (Heteroscedastic)');
xlim([noises(1), noises(end)])