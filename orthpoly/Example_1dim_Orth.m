% Initializing Random Number Generator (octave support)
rng_i shuffle;
close all;

% function to be reconstructed
func = @(x)4 - 0.6*x - x.^2 + 1.7*x.^3;

disp('Reconstructing function'); disp(func);
a = -1.6; b = 2.1; % interval
h = 0.02; % step

% initial function
x = a:h:b;
y = func(x);

% taking N data points
N = 5;
rx = affine_transform(rand(1, N), [a; b], [0; 1]);
ry = func(rx);

% plotting points
figure(1);
subplot(211);
plot(x, y, '-.g', rx, ry, 'ok');
hold on; grid on;

% constructing orthogonal polynomials
deg = 3;
c01 = [0; 1];
F = orthpoly(deg, 1, c01(1), c01(2));
sigma = deglexord(deg, 1);
E = eye(deg + 1);

% orthogonality test - result F must be close to identity
orthogonality_test(F, deg, 1, c01(1), c01(2), 1e-12)

% getting values
R = EvalPolyOrth(E, rx', sigma, F, sigma);

coefs = (R'*R)\R'*ry';
disp('Coefficients are'); disp(coefs);
eq = disp(coefs');

% finding values of a reconstructed function
y1 = (EvalPolyOrth(E, x', sigma, F, sigma) * coefs)';

disp('Reconstructed function error:');
err = norm(y - y1)

% plotting reconstructed function
plot(x, y1, 'm');
legend('Initial function', 'Data points', 'Reconstructed function');
title(sprintf('Function reconstruction with %u points (order %u)\nCoefficients are %s', N, deg, eq));
clear eq;

% error
subplot(212);
plot(x, abs(y - y1), 'r');
title('Error');
