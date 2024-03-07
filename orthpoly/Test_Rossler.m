rng_i shuffle;
close all;
warning off;

% Rossler system Simulation
Tmax = 50; % Time end
h = 0.01; % Step
start_point = [4 -2 0]; % Initial point
disp('Start point is '), disp(start_point)

[t, x] = ode45(@Rossler, 0:h:Tmax, start_point);
y = transpose(Rossler(0, x'));

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

% Rossler system has fc = 3 variables (x, y, z) and every equation polynom degree is less or equal to 2
fc = size(y, 2);
deg = 2;
sigma = berndeg(deg, fc);

% Applying affine transform to fit the [0, 1] interval
c01 = repmat([0; 1], 1, fc);
[tx, c] = affine_transform(rx, c01);

% Finding Bernstein monomial values in every point
B = bernbase(tx, sigma);
coefs = (B'*B)\B'*ry; % LSM solution

disp('Bernstein coefficients are:'); disp(coefs);
disp('Domain:'); disp(c);

% Testing system on real values
tx1 = affine_transform(x, c01, c);
y1 = bernbase(tx1, sigma) * coefs;

disp('Reconstructed system error:');
err = norm(y - y1)

% Solving scary ODE in Bernstein polynomials
[~, x1] = ode45(@(t,x)transpose(bernbase(affine_transform(x', c01, c), sigma) * coefs), 0:h:Tmax, start_point);

plot3(x1(:,1), x1(:,2), x1(:,3), '-y', 'DisplayName', 'Reconstructed function');
title(sprintf('Rossler system reconstruction\nError: %g', err));
