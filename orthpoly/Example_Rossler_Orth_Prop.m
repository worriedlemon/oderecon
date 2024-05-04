close all;
warning off;

% System Simulation
start_point = [4 -2 0]; % Initial point
Tmax = 50; % Time nd
h = 0.01; % Step

[t, x] = ode45(@Rossler, 0:h:Tmax, start_point);

deg = 2; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = size(y, 2); % Equations count
sigma = deglexord(deg, vc);

F = orthpoly_t(sigma, t, x); % Getting orthogonal polynomials matrix
mc = size(F, 1); % Monomials count

coefs = zeros(mc, eqc);
for j = 1:mc
    E = EvalPoly(F(j, :)', x, sigma);
    for i = 1:eqc
        coefs(j, i) = trapz(x(:, i), E);
    end
end


disp("\nCoefficients:");
coefs = F' * coefs %#ok

H = mat2cell(coefs, mc, ones(1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

figure(1);
plot3(x(:, 1), x(:, 2), x(:, 3), 'b');
hold on; grid on;

[~, x1] = ode45(@(t, x)oderecon(H, T, t, x), t, start_point);
plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'r');