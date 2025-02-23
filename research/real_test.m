close all;
warning off;

% System Simulation
sysname = 'SprottB';
sys_f = 'P7_1.txt';

lsm_on = 1;
phase_on = 1;
sync_on = 1;

disp(['Reconstructing ', sysname, ' system, Variant ', sys_f]);
x_r = readmatrix(['research/IGN_', sysname, '/', sys_f]);
h_pred = 6e-6;

h = h_pred; 

t = h * x_r(:, 1);
x_r = x_r(:, 2 * (1:3)) ./ [1 (7.07 * 10000 / 120) 1];

start_point = x_r(1, :); % Initial point
[~, x_orig] = ode45(@SprottB, t, start_point);

deg = 3; % Degree of reconstructed function
vc = size(x_r, 2); % Variables count
eqc = vc; % Equations count
sigma = deglexord(deg, vc);
mc = size(sigma, 1); % Monomials count

F = orthpoly_t(sigma, t, x_r) %#ok Getting orthogonal polynomials matrix

x = sgolayfilt(x_r, deg, 2 * deg + 1);

coefs = zeros(mc, eqc);
coefs_reg = coefs;
% E = EvalPoly(F', x, sigma);
% for j = 1:mc
%     for i = 1:eqc
%         coefs(j, i) = trapz(x(:, i), E(:, j));
%     end
% end

y = diff4(x, t);
tol = 1e-1;
for j = 1:vc
    [coefs(:, j), ~, ~, coefs_reg(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
end
coefs
disp('Regular Coefficients:'); Ho = coefs_reg %#ok

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
    plot3(x(:, 1), x(:, 2), x(:, 3), 'b', 'DisplayName', 'Smoothed original data');
    hold on; grid on;
    %plot3(x_orig(:, 1), x_orig(:, 2), x_orig(:, 3), 'g', 'DisplayName', 'Simulated data');
    %plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'r', 'DisplayName', 'Orthogonal Polynomials reconstruction');
    if (lsm_on)
        %plot3(x2(:, 1), x2(:, 2), x2(:, 3), 'g', 'DisplayName', 'LSM');
        plot3(x2(:, 1), x2(:, 2), x2(:, 3), 'r', 'DisplayName', 'LSM');
    end
    xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$')
    set(gca,'TickLabelInterpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    legend show
end

if (sync_on)
    sync_coef = 1e5;
    xo_slave = RK4SyncFOH(x, Ho_c, T, t, sync_coef);
    
    figure(2);
    semilogy(t, vecnorm(x - xo_slave, 2, 2), 'Color', [0 0 1 0.2]);
    hold on; grid on;
    semilogy(t, (t * 0) + mean(vecnorm(x - xo_slave, 2, 2)), 'Color', [0 0 1 1], 'LineWidth', 3);
    legend('', 'Orthpoly (Mean)');
    if (lsm_on)
        xt_slave = RK4SyncFOH(x, Ht_c, T, t, sync_coef);
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