close all;
warning off;

% System Simulation
sys = @FHN;

% Use delMinorTerms
delminor = 0;

start_point = [0.1 0.1]; % Initial point
Tmax = 25; % Time end
h = 1e-2; % Step

[t, x] = ode45(sys, 0:h:Tmax, start_point);

deg = 3; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = vc; % Equations count
sigma = deglexord(deg, vc);

F = orthpoly_t(sigma, t, x); % Getting orthogonal polynomials matrix
mc = size(F, 1); % Monomials count

coefs = zeros(mc, eqc);
if ~delminor
    E = EvalPoly(F', x, sigma);
    
    for j = 1:mc
        for i = 1:eqc
            coefs(j, i) = trapz(x(:, i), E(:, j));
        end
    end
else
    % We need to know derivatives
    y = [diff(x) / h; (x(end, :) - x(end - 1, :)) / h]; % first order
    
    for j = 1:vc
        [coefs(:, j), ~] = delMinorTerms_dx(x(:, j), x, y(:, j), F, sigma, 2e-3, 0);
    end
end

disp('Orthogonal Coefficients:'); coefs %#ok
disp('Regular Coefficients:'); coefs = F' * coefs %#ok

H = mat2cell(coefs, mc, ones(1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

figure(1);
plot(x(:, 1), x(:, 2));
hold on; grid on;

[~, x1] = ode45(@(t, x)oderecon(H, T, t, x), t, start_point);
plot(x1(:, 1), x1(:, 2), 'r');
xtickformat('$%g$'); ytickformat('$%g$');
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
%title('FitzHugh-Nagumo model reconstruction');
legend('Initial system', 'Reconstructed system')