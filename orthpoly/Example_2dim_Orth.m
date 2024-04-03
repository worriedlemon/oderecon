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
[x,y] = meshgrid(a:h:b, a:h:b);
f = func(x,y);

figure(1);
subplot(211);
surf(x, y, f);
colormap([0.6, 1, 1]);
hold on; grid on;

% taking N data points
N = 40;
deg = 4; % polynom degree

rx = affine_transform(rand(1, N), [a; b], [0; 1]);
ry = affine_transform(rand(1, N), [a; b], [0; 1]);
rf = func(rx, ry);

% plotting points
scatter3(rx, ry, rf, 35, 'marker', 'o', 'markeredgecolor', 'red', 'markerfacecolor', 'white', 'linewidth', 2);
legend('Initial function', 'Data points')
title(sprintf('Function reconstruction with %u points (order %u)', N, deg));

% constructing orthogonal basis
deg = 4;
vc = 2;
c01 = [0; 1];
F = orthpoly(deg, vc, c01(1), c01(2));
sigma = deglexord(deg, vc);
E = eye(size(sigma, 1));

% Orthogonality test -- must be close to identity
orthogonality_test(F, deg, vc, c01(1), c01(2));

% getting values
R = EvalPolyOrth(E, [rx; ry]', F, sigma);
coefs = (R'*R)\R'*rf';
disp('Coefficients are'); disp(coefs);

% finding values of a reconstructed function
f1 = zeros(size(x));
for i = 1:size(x,1)
    f1(i,:) = EvalPolyOrth(E, [x(i,:); y(i,:)]', F, sigma) * coefs;
end

% finding error
disp('Reconstructed function error:')
err = norm(f - f1)

% plotting reconstructed function
subplot(212);
surf(x, y, f1);
hold on; grid on;
scatter3(rx, ry, rf, 35, 'marker', 'o', 'markeredgecolor', 'red', 'markerfacecolor', 'white', 'linewidth', 2);
legend('Reconstructed function', 'Data points');
title(sprintf('Error: %g', err))
