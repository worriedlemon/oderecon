rng_i default;
close all;
warning off;

sys = @Lorenz;
func = @(x,y)(x.^2 + y.^2 - 11).^2 + (x + y.^2 - 7).^2;

start_point = [4 -2 0]; % Initial point
Tmax = 50; % Time end
h = 0.01; % Step

%[~, x] = ode45(sys, 0:h:Tmax, start_point);
%y = transpose(sys(0, x'));
a = -4; b = 4;
span = (a:h:b)';
[x1, x2] = meshgrid(span);
x = [reshape(x1, length(span)^2, 1), reshape(x2, length(span)^2, 1)];
clear x1 x2;
y = func(x(:, 1), x(:, 2));

deg = 4; % Degree of reconstructed function
vc = size(x, 2); % Variables count
eqc = size(y, 2); % Equations count

% Taking N random points
N = 10000;
ri = randperm(size(y, 1), N);
rx = x(ri,:);
ry = y(ri,:);

% Normalizing values to fit [0; 1]
c01 = repmat([0; 1], 1, vc);
[tx, cfrom] = affine_transform(rx, c01);

% Dividing space on cubes
Nd = 50;
delta = [-1, 1] * c01 / Nd; % Deltas for cube
J = prod(delta / 2, 2); % Jacobian for quadrature

sigma = deglexord(deg, vc); % Ordered monomials
linapprox = sigma(1:vc + 1, :); % Monomials for linear approximation

mc = size(sigma, 1); % Monomials count
[F, norms] = orthpoly(deg, vc, 0, 1, 1e-9, 0); % Getting orthogonal polynomials matrix and norms
opt = {'orth', F, sigma};

coefs = zeros(mc, eqs);
for eq = 1:eqc
    for i = 1:mc % For each monomial
        for k = 0:Nd^vc - 1 % For each elementary boundary
            
            % Workaround for nested cycles to be as the single one
            j = dec2base_imp(k + Nd^vc, Nd);
            j = j(2:vc + 1);
            
            % Taking next bounded cube
            bounds = [j .* delta; (j + 1) .* delta];
            
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
coefs = coefs ./ (norms .^ 2); % Coefficients for orthogonal monomials

H = mat2cell(coefs, [mc], repmat([1], 1, eqc));
T = mat2cell(repmat(sigma, 1, eqs), [mc], repmat([vc], 1, eqc));

span = (a:0.1:b)';
[x1, x2] = meshgrid(span);
x = [reshape(x1, length(span)^2, 1), reshape(x2, length(span)^2, 1)];
clear x1 x2;
y = func(x(:, 1), x(:, 2));

%[xc, cfrom] = affine_transform(x, c01, cfrom);
%start_point_c = affine_transform(start_point, c01, cfrom);
%[~, xc1] = ode45(@(t, x)oderecon(H, T, t, affine_transform(x', c01)', opt), 0:h:Tmax, start_point_c);
y1 = EvalPoly(coefs, affine_transform(x, c01, cfrom), sigma, opt);

figure(1);
plot3(x(:, 1), x(:, 2), y, '.b');
zlim([-50; 900]);

figure(2);
plot3(x(:, 1), x(:, 2), y1, '.r');
zlim([-50; 900]);

prettyOrth(H, T, F, sigma, 0);