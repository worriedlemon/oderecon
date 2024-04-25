rng_i default;
close all;
warning off;

sys = @(x,y)(x.^2 + y.^2 - 11).^2 + (x + y.^2 - 7).^2;

h = 0.1; % Step
a = -4;
b = 4;
span = (a:h:b)';
[x1, x2] = meshgrid(span);
x = [reshape(x1, length(span)^2, 1), reshape(x2, length(span)^2, 1)];
clear x1, x2;
y = sys(x(:, 1), x(:, 2));

deg = 4; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = size(y, 2); % Equations count

% Taking N random points
N = 5000;
ri = randperm(size(y, 1), N);
rx = x(ri,:);
ry = y(ri,:);

% Normalizing values to fit [0; 1]
c01 = repmat([0; 1], 1, vc);
[tx, cfrom] = affine_transform(rx, c01);

% Dividing space on cubes
Nd = 32;
delta = [-1, 1] * c01 / Nd; % Deltas for cube
J = prod(delta / 2, 2); % Jacobian for quadrature

sigma = deglexord(deg, vc); % Ordered monomials
linapprox = sigma(1:vc + 1, :); % Monomials for linear approximation

mc = size(sigma, 1); % Monomials count
[F, norms] = orthpoly(deg, vc, 0, 1, 1e-9, 0); % Getting orthogonal polynomials matrix and norms
opt = {'orth', F, sigma};

coefs = zeros(mc, eqc);
for eq = 1:eqc
    for i = 1:mc % For each monomial
        for k = 0:Nd^vc - 1 % For each elementary boundary
            
            % Workaround for nested cycles to be as the single one
            if (Nd < 2)
                j = zeros(1, vc);
            else
                j = dec2base_imp(k + Nd^vc, Nd);
                j = j(2:vc + 1);
            end
            
            % Taking next bounded cube
            bounds = [j; j + 1] .* delta;
            
            % Finding indexes of points which fall inside of bounded cube
            idx = find(prod(tx <= bounds(2, :) & tx >= bounds(1, :), 2));
            
            if length(idx) == 0
                w = 1;
                E = 0;
            else
                if length(idx) < vc + 1
                    % Rectangle mathod
                    w = 2^vc;
                    E = mean(sum(ry(idx, eq) .* EvalPoly(F(i, :)', tx(idx, :), sigma), 2)); % Average value of f * g_i
                else
                    % Trapezoid method
                    w = 1;
                    B = EvalPoly(eye(vc + 1), tx(idx, :), linapprox);
                    pl = (B'*B)\B'*ry(idx, eq); % LSM for finding approximation plane
                    ps = bounds(1, :) + (dec2bin(transpose(0:2^vc - 1)) - '0') .* delta; % Magic for finding vertex points
                    E = sum(EvalPoly(pl, ps, linapprox) .* EvalPoly(F(i, :)', ps, sigma)); % Sum of values f * g_i
                end
            end
            coefs(i, eq) = coefs(i, eq) + w * J * E; % Generalized Riemann sum
        end
    end
end
coefs = coefs ./ (norms .^ 2); % Coefficients for orthogonal monomials (normalized function f(t), t = (x - a)/(b - a))

H = mat2cell(coefs, [mc], repmat([1], 1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), [mc], repmat([vc], 1, eqc));

x1 = affine_transform(x, c01, cfrom);
y1 = EvalPoly(F' * coefs, x1, sigma);

figure(1);
plot3(x(:, 1), x(:, 2), y, '.b');

figure(2);
plot3(x1(:, 1), x1(:, 2), y1, '.r');

prettyOrth(H, T, F, sigma, 0);