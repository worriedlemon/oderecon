close all;
warning off;

% System Simulation
sysname = 'SprottB';
sys_f = 'P8_1.txt';

lsm_on = 0;
phase_on = 1;
sync_on = 0;

disp(['Reconstructing ', sysname, ' system, Variant ', sys_f]);
x_r = readmatrix(['research/IGN_', sysname, '/', sys_f]);
h_pred = 1; % not a real step, for compatibility

h = h_pred; 

t = h * x_r(:, 1);
x_r = x_r(:, 2 * (1:3));

start_point = x(1, :); % Initial point
Tmax = 100; % Time end

deg = 3; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = vc; % Equations count
sigma = deglexord(deg, vc);

F = orthpoly_t(sigma, t, x) %#ok Getting orthogonal polynomials matrix
mc = size(F, 1); % Monomials count

x = sgolayfilt(x_r, deg, 2 * deg + 1);

coefs = zeros(mc, eqc);
E = EvalPoly(F', x_r, sigma);
for j = 1:mc
    for i = 1:eqc
        coefs(j, i) = trapz(x(:, i), E(:, j));
    end
end

disp('Regular Coefficients:'); Ho = F' * coefs %#ok

Ho_c = mat2cell(Ho, mc, ones(1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

[~, x1] = ode45(@(t, x)oderecon(Ho_c, T, t, x), t, start_point);

if (lsm_on)
    y = [diff(x); (x(end, :) - x(end - 1, :))] / h; % first order
    B = EvalPoly(eye(mc), x, sigma);
    delta = 0.01;
    
    disp('Regular Coefficients (LSM):');
    Ht = (B'*B + delta*eye(mc))\B'*y %#ok
    Ht_c = mat2cell(Ht, mc, ones(1, eqc));

    [~, x2] = ode45(@(t, x)oderecon(Ht_c, T, t, x), t, start_point);
end

if (phase_on)
    figure(1);
    plot3(x_r(:, 1), x_r(:, 2), x_r(:, 3), 'g--', 'DisplayName', 'Original data');
    hold on; grid on;
    plot3(x(:, 1), x(:, 2), x(:, 3), 'b', 'DisplayName', 'Smoothed original data');
    plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'r', 'DisplayName', 'Orthogonal Polynomials');
    if (lsm_on)
        plot3(x2(:, 1), x2(:, 2), x2(:, 3), 'g', 'DisplayName', 'LSM');
    end
    xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$')
    set(gca,'TickLabelInterpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    legend show
end

if (sync_on)
    N = length(t);
    xo_slave = RK4SyncFOH(x, Ho_c, T, t);
    
    figure(2);
    semilogy(t, vecnorm(x - xo_slave, 2, 2), 'Color', [0 0 1 0.2]);
    hold on; grid on;
    semilogy(t, (t * 0) + mean(vecnorm(x - xo_slave, 2, 2)), 'Color', [0 0 1 1], 'LineWidth', 3);
    legend('', 'Orthpoly (Mean)');
    if (lsm_on)
        xt_slave = RK4SyncFOH(x, Ht_c, T, t);
        semilogy(t, vecnorm(x - xt_slave, 2, 2), 'Color', [1 0 0 0.2]);
        semilogy(t, (t * 0) + mean(vecnorm(x - xt_slave, 2, 2)), 'Color', [1 0 0 1], 'LineWidth', 3);
        legend('', 'Orthpoly (Mean)', '', 'LSM (Mean)');
    end
    legend show;
    xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
    set(gca, 'TickLabelInterpreter', 'latex');
    %xlabel('Time $t$, s', 'Interpreter', 'latex');
    ylabel('Synchronization error $\overline{\zeta}$', 'Interpreter', 'latex');
    xlim([t(end - 1000), t(end)]);
end