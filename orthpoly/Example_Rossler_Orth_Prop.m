rng_i default;
close all;
warning off;

% Rossler system Simulation
start_point = [4 -2 0]; % Initial point
Tmax = 10; % Time nd
h = 0.1; % Step

[t, x] = ode45(@Rossler, 0:h:Tmax, start_point);
y = transpose(Rossler(0, x'));

deg = 2; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = size(y, 2); % Equations count

[F, nrms] = orthpoly_t(deg, vc, t, x, 0); % Getting orthogonal polynomials matrix and norms
mc = size(F, 1); % Monomials count

coefs = zeros(mc, eqc);
for j = 1:mc
    coefs(j, :) = trapz(t, EvalPoly(transpose(F(j, :) / nrms(j)), x, sigma) .* y);
end

coefs = (F ./ nrms)' * coefs;

disp("\nCoefficients:");
coefs

H = mat2cell(coefs, [mc], repmat([1], 1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), [mc], repmat([vc], 1, eqc));

figure(1);
plot3(x(:, 1), x(:, 2), x(:, 3), 'b');
hold on; grid on;

[~, x1] = ode45(@(t, x)oderecon(H, T, t, x), t, start_point);
plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'r');