% Initializing Random Number Generator (octave support)
rng_i default;
close all;

% function to be reconstructed (Himmelblau)
%func = @(x,y)170 - 14*x - 21*x.^2 - 36*y.^2 + 2*x.*y.^2 + x.^4 + 2*x.^2.*y.^2 + 2*y.^4;
func = @(x,y)(x.^2 + y.^2 - 11).^2 + (x + y.^2 - 7).^2;

disp('Reconstructing function'); disp(func);
a = -4; b = 4; % interval
h = 0.1; % step

% initial function plot
x = a:h:b;
y = x;
f = zeros(length(x), length(y));
for i = 1:length(x)
  for j = 1:length(y)
    f(i, j) = func(x(i), y(j));
  end
end

figure(1);
subplot(211);
plot3(x, y, f, '.g');
hold on; grid on;

% taking N data points
N = 40;
deg = 4; % polynom degree

rx = rand(1, N) * (b - a) + a;
ry = rand(1, N) * (b - a) + a;
rf = func(rx, ry);

% plotting points
scatter3(rx, ry, rf, 'k');
legend('Initial function', 'Data points')
title(sprintf('Function reconstruction with %u points (order %u)', N, deg));

% finding Bernstein base monomials values
c01 = [0 0; 1 1];
[txy, c] = affine_transform([rx; ry]', c01);
deg = 4;

sigma = berndeg(deg, 2);
E = eye(size(sigma, 1));
B = EvalPoly(E, txy, sigma, 'bernstein');

coefs = (B'*B)\B'*rf';
disp('Bernstein coefficients are'); disp(coefs);

disp('Domain:'); disp(c);

% finding values of a reconstructed function
f1 = zeros(length(x), length(y));
for i = 1:length(x)
  for j = 1:length(y)
    txy1 = affine_transform([x(i); y(j)]', c01, c);
    f1(i, j) = EvalPoly(E, txy1, sigma, 'bernstein') * coefs;
  end
end

% finding error
disp('Reconstructed function error:')
err = norm(f - f1)

% plotting reconstructed function
subplot(212);
plot3(x, y, f1, '.m');
hold on; grid on;
scatter3(rx, ry, rf, 'k');
legend('Reconstructed function', 'Data points');
title(sprintf('Error: %g', err))
