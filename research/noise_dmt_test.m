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
start_point = x(end, :);
%then, go experiment
[t, x] = ode78(system, 0:h:Tmax, start_point);


y = transpose(system(0, x'));
[N, vc] = size(x);

eqc = size(y, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

noises = logspace(-4, 1, 10);
tols = logspace(-3, 0, 20);
erro = zeros(length(noises), length(tols));
mcs = erro;
for k = 1:length(noises)
    noise_amp = noises(k);
    rx = x + noise_amp * randn(N, vc);
    ry = diff4(rx, t); % fourth order
    
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    for i = 1:length(tols)
        tol = tols(i);
        for j = 1:eqc
            [~, ~, ~, Ho(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
        end
        mcs(k, i) = sum(Ho ~= 0, "all");
        erro(k, i) = norm(Ho - Href);
    end
end

figure(1);
cls = colormap("hsv");
cls = cls(floor(linspace(1, size(cls, 1), 6/5 * length(noises))), :);
subplot(2, 1, 1);
for k = 1:length(noises)
    nstr = num2str(noises(k));
    loglog(tols, erro(k, :), Marker='.', Color=cls(k, :));
    hold on; grid on;
end
xtickformat('$%g$'); ytickformat('$%g$');
xlabel('Tolerance $\eta$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

subplot(2, 1, 2);
for k = 1:length(noises)
    nstr = num2str(noises(k));
    semilogx(tols, mcs(k, :), Marker='.', Color=cls(k, :), DisplayName=['$\sigma=', nstr, '$']);
    hold on; grid on;
end
loglog(tols, tols * 0 + sum(Href ~= 0, "all"), 'k--', DisplayName='Actual monomials count');
legend(Interpreter='latex');
xtickformat('$%g$'); ytickformat('$%g$');
xlabel('Tolerance $\eta$', 'Interpreter', 'latex');
ylabel('Monomials count $L$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ylim([0; mc * vc]);

sgtitle('Monomials count depending on tolerance with different levels of noise magnitude $\sigma$', Interpreter='Latex')