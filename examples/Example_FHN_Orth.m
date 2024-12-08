close all;
warning off;

% System Simulation
sys = @FHN;

% Use delMinorTerms
delminor = 1;

start_point = [0.1 0.1]; % Initial point
Tmax = 25; % Time end
h = 1e-2; % Step

[t, x] = ode45(sys, 0:h:Tmax, start_point);

deg = 3; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = vc; % Equations count
sigma = deglexord(deg, vc);

F = orthpoly_t(sigma, t, x) %#ok Getting orthogonal polynomials matrix
mc = size(F, 1); % Monomials count

coefs = zeros(mc, eqc);
coefs_reg = zeros(mc, eqc);

if ~delminor
    E = EvalPoly(F', x, sigma);
    for j = 1:mc
        for i = 1:eqc
            coefs(j, i) = trapz(x(:, i), E(:, j));
        end
    end
    coefs_reg =  F' * coefs;
else
    % We need to know derivatives
    y = diff4(x,t);
    
    tol = 1e-3; % tolerance
    for j = 1:vc
        [coefs(:, j), ~, ~, coefs_reg(:, j)] = delMinorTerms_dy(t,x(:, j), x, y(:, j), F, sigma, tol, 0);
    end
end


disp('Orthogonal Coefficients:'); coefs_orth = coefs %#ok
disp('Regular Coefficients:'); coefs = coefs_reg %#ok

T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));
H = mat2cell(coefs, mc, ones(1, eqc));


% Display equations in orthogonal polynomials
H_orth = mat2cell(coefs_orth, mc, ones(1, eqc));
incf = 1;
prettyOrth(H_orth,T,F,sigma,incf);

figure(1);
plot(x(:, 1), x(:, 2));
hold on; grid on;

[~, x1] = ode45(@(t, x)oderecon(H, T, t, x), t, start_point);
plot(x1(:, 1), x1(:, 2), 'r');
xtickformat('$%g$'); ytickformat('$%g$');
set(gca,'TickLabelInterpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
legend('Initial system', 'Reconstructed system');

prettyABM(H,T)
