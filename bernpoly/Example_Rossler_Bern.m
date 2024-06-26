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

% Rossler system has fc = 3 variables (x, y, z) and every equation polynom degree is less or equal to 2
fc = size(y, 2);
deg = 2;

% Taking N data points
N = 30;
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

% Getting order ideal
sigma = berndeg(deg, fc);

% Applying affine transform to fit the [0, 1] interval
c01 = repmat([0; 1], 1, fc);
[tx, c] = affine_transform(rx, c01);

% Use LSM for fitting the equations with the proper coefficients
eta = 1e-5;
H = cell(1, fc);
T = cell(1, fc);

% Reconstruct each equation
ry0 = ry;
for i = 1:3
    [hi, tau] = delMinorTerms(tx, ry(:,i), sigma, eta, 'bernstein'); % Get equation and basis
    ry0(:,i) = EvalPolyBern(hi, tx, tau); % Get values
    
    H{1,i} = hi;
    T{1,i} = tau;
end

disp('Domain:'); disp(c);
disp('System reconstruction error:');
err = vecnorm(ry - ry0)

% Solving scary ODE in Bernstein polynomials
[~, x1] = ode45(@(t,x)oderecon(H, T, t, affine_transform(x', c01, c)', 'bernstein'), 0:h:Tmax, start_point); %solve ODE

plot3(x1(:,1), x1(:,2), x1(:,3), '-y', 'DisplayName', 'Reconstructed function');
title(sprintf('Rossler system reconstruction\nError: %s', disp(err)));

% Output as equations, saving equations for future processing
equations = prettyBernstein(H, T, 0)