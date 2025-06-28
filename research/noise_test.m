rng_i default;
warning off;
close all

% System coefficients
Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.3 0 0 0 0 0 0 0; 0.3 0 0 -5.7 0 0 1 0 0 0]';

% Used system
%system = @Rossler
system = @Lorenz

Href = eval(['H', lower(func2str(system))]); % Used coefficients
Ttrans = 100; % Transient time
Tmax = 5; % Time of experimental series for Lorenz
%Tmax = 10; % Time of experimental series for Rossler
h = 1e-2; % Step
start_point = [4 -2 0]; % Initial point

% %first, go transient
[~, x] = ode78(system, 0:h:Ttrans, start_point);
start_point = x(end,:);
%then, go experiment
[t, x] = ode78(system, 0:h:Tmax, start_point);

y = transpose(system(0, x'));

[N, vc] = size(x);

eqc = size(y, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

delta = 1e-7;  % regularization parameter
tol = 5e-2;    % DMT tolerance
lambda = 1e-1; % SINDy sparsification parameter

noises = logspace(-5, 1, 200);

% Homoscedasticity test
errt = []; erro = []; errs = []; erro2 = [];
for noise_amp = noises
    rx = x + noise_amp * randn(N, vc);
    ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
    %ry = diff2(rx) / h; % second order
    %ry = diff4(rx, t); % fourth order
    
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc
        for j = 1:mc
            %Ho(j, i) = trapz(t, E(:, j) .* ry(:, i));
            %Ho(j, i) = trapz(rx(:, i), E(:, j));
            %Ho(j, i) = E(end, j)*rx(end, i) - E(1, j)*rx(1, i) -  integrate_simpvar(E(:, j), rx(:, i));


            %Ho(j, i) = integrate_simpvar(t, E(:, j) .* ry(:, i));
            Ho(j, i) = intdiff4(rx(:, i),E(:, j));
           
            
            %Ho(j, i) = integrate_simpvar(rx(:, i), E(:, j));
            % tmp = integrate_bool(E(:, j) .* ry(:, i),0,h);
            % Ho(j, i) = tmp(end);
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B + delta*eye(mc))\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];

    Ht1dmt = Ht1;
    tol = 5e-2;
    for j = 1:vc
        [~, ~, ~, Ht1dmt(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
    end
    erro2 = [erro2 norm(Ht1dmt - Href)];
    
    Ht2 = sindy(rx, ry, sigma, lambda);
    errs = [errs norm(Ht2 - Href)];
end

figure(1);
loglog(noises, errt, 'r', DisplayName='LSM (Homoscedastic)');
hold on; grid on;
loglog(noises, errs, 'g', DisplayName='SINDy (Homoscedastic)');
loglog(noises, erro, 'b', DisplayName='OrthPoly (Homoscedastic)');
loglog(noises, erro2, 'm', DisplayName='OrthPoly DMT (Homoscedastic)');
%title(['Noise resistance (', func2str(system), ')']);
legend show;
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlim([noises(1), noises(end)])

% Heteroscedasticity test
errt = []; erro = []; errs = []; erro2 = [];
for noise_amp = noises
    dmx = x - mean(x);
    rx = x + noise_amp * randn(N, vc) .* dmx;
    ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h];
    %ry = diff2(rx) / h; % second order
    %ry = diff4(rx, t); % fourth order
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc 
        for j = 1:mc
            %Ho(j, i) = trapz(rx(:, i), E(:, j));
            %Ho(j, i) = integrate_simpvar(t, E(:, j) .* ry(:, i));
            %Ho(j, i) = intdiff2(rx(:, i),E(:, j));
            %Ho(j, i) = trapz(t, E(:, j) .* ry(:, i));
            Ho(j, i) = intdiff4(rx(:, i),E(:, j));
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B)\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];

    Ht1dmt = Ht1;
    tol = 5e-2;
    for j = 1:vc
        [~, ~, ~, Ht1dmt(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
    end
    erro2 = [erro2 norm(Ht1dmt - Href)];
    
    Ht2 = sindy(rx, ry, sigma, lambda);
    errs = [errs norm(Ht2 - Href)];
end

figure(2);
loglog(noises, errt, 'r', DisplayName='LSM (Heteroscedastic)');
hold on; grid on;
loglog(noises, errs, 'g', DisplayName='SINDy (Heteroscedastic)');
loglog(noises, erro, 'b', DisplayName='OrthPoly (Heteroscedastic)');
loglog(noises, erro2, 'm', DisplayName='OrthPoly DMT (Heteroscedastic)');
%title(['Noise resistance (', func2str(system), ')']);
legend show;
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlim([noises(1), noises(end)])