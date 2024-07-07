close all;
warning off;

% System Simulation
sys = @Lorenz;
sysname = 'Lorenz';
start_point = [4 -2 0]; % Initial point
Tmax = 50; % Time nd
h = 1e-2; % Step

[t, x] = ode45(sys, 0:h:Tmax, start_point);

deg = 2; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = vc; % Equations count
sigma = deglexord(deg, vc);

F = orthpoly_t(sigma, t, x); % Getting orthogonal polynomials matrix
mc = size(F, 1); % Monomials count

E = EvalPoly(F', x, sigma);
coefs = zeros(mc, eqc);
for j = 1:mc
    for i = 1:eqc
        coefs(j, i) = trapz(x(:, i), E(:, j));
    end
end


disp("\nCoefficients:");
coefs = F' * coefs %#ok

H = mat2cell(coefs, mc, ones(1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

figure(1);
plot3(x(:, 1), x(:, 2), x(:, 3));
hold on; grid on;

[~, x1] = ode45(@(t, x)oderecon(H, T, t, x), t, start_point);
plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'r');
xlabel('\itx'); ylabel('\ity'); zlabel('\itz');
title([sysname, ' dynamical system reconstruction']);

legend('Initial system', 'Reconstructed system')