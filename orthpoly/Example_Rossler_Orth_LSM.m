rng_i default;
close all;
warning off;

% Rossler system Simulation
Tmax = 50; % Time end
h = 0.01; % Step
start_point = [4 -2 0]; % Initial point
disp('Start point is '), disp(start_point)

[t, x] = ode45(@Rossler, 0:h:Tmax, start_point);
y = transpose(Rossler(0, x'));

% Rossler system has vc = 3 variables (x, y, z) and every equation polynom degree is less or equal to 2
vc = size(y, 2);
deg = 2;

% Taking N data points
N = 15;
ri = randi(size(y, 1), 1, N);
rx = x(ri,:);
ry = y(ri,:);

% Plotting initial Rossler system and data points
figure(1);
plot3(x(:,1), x(:,2), x(:,3), 'b', start_point(1), start_point(2), start_point(3), 'r*');
hold on;
scatter3(rx(:,1), rx(:,2), rx(:,3), 'k');
xlabel('\itx'); ylabel('\ity'); zlabel('\itz');
legend('Rossler attractor', 'Start point', 'Random points');

% Getting order ideal, constructing orthogonal polynomials basis
sigma = deglexord(deg, vc);
[F, nrms] = orthpoly(deg, vc, 0, 1, 1e-9, 0);
orthogonality_test(F, deg, vc, 0, 1, 'brief', 1e-6, nrms);

opt = {'orth', F, sigma};

% Use LSM for fitting the equations with the proper coefficients
eta = 1e-7;
H = cell(1, vc);
T = cell(1, vc);

% Reconstruct each equation
ry0 = ry;
E = EvalPolyOrth(eye(10), rx, sigma, F, sigma);
h0 = (E'*E)\E'*ry;
for i = 1:3
    [hi, tau] = delMinorTerms(rx, ry(:,i), sigma, eta, opt, h0(:, i)); % Get equation and basis
    ry0(:,i) = EvalPolyOrth(hi, rx, tau, F, sigma); % Get values
    
    H{1,i} = hi;
    T{1,i} = tau;
end

disp('System reconstruction error:');
err = vecnorm(ry - ry0)

% Solving ODE
[~, x1] = ode45(@(t,x)oderecon(H, T, t, x, opt), 0:h:Tmax, start_point);

plot3(x1(:,1), x1(:,2), x1(:,3), '-y', 'DisplayName', 'Reconstructed function');
title(sprintf('Rossler system reconstruction\nError: %s', disp(err)));

% Output as equations, saving equations for future processing
equations = prettyOrth(H, T, F, sigma, 0)