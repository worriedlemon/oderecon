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

% finding Bernstein base monomials values
eps = 0.1;
[tx, c] = affine_transform(rx, [0; 1]);
deg = 3;

sigma = berndeg(deg, 1);
B = bernbase(tx', sigma);

coefs = (B'*B)\B'*ry';
disp('Bernstein coefficients are'); disp(coefs);
eq = disp(coefs');

disp('Domain:'); disp(c');

% finding values of a reconstructed function
y1 = (bernbase(affine_transform(x, [0; 1], c)', sigma) * coefs)';

disp('Reconstructed function error:');
err = norm(y - y1)

% plotting reconstructed function
plot(x, y1, 'm');
legend('Initial function', 'Data points', 'Reconstructed function');
title(sprintf('Function reconstruction with %u points (order %u)\nBernstein coefficients are %s', N, deg, eq));
clear eq;

% error
subplot(212);
plot(x, abs(y - y1), 'r');
title('Error');
