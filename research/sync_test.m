rng_i default;
warning off;

Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.2 0 0 0 0 0 0 0; 0.2 0 0 -5.7 0 0 1 0 0 0]';
system = @Lorenz;
sysname = 'Lorenz';
Href = Hlorenz;

Tmax = 100; % Time end
h = 5e-3; % Step
start_point = [0.1 0 0.1]; % Initial point

[t, x] = ode45(system, 0:h:Tmax, start_point);

[N, vc] = size(x);
eqc = vc;
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

rx = x;
ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
%ry = diff2(rx) / h; % second order
%ry = diff4(rx, t); % fourth order


B = EvalPoly(eye(mc), rx, sigma);
Ht = (B'*B)\B'*ry;

F = orthpoly_t(sigma, t, rx);

Ho = zeros(mc, eqc);
E = EvalPoly(F', rx, sigma);
for i = 1:eqc
    for j = 1:mc
        Ho(j, i) = trapz(rx(:, i), E(:, j));
    end
end
Ho = F' * Ho;

Tmax = 100;
h = 1e-3;
t = 0:h:Tmax;
N = length(t);

[~, x] = ode45(system, t, start_point);
xt_slave = [start_point; zeros(N - 1, vc)];
xo_slave = xt_slave;
sync_coef = 2 * [1 1 1];
for i = 1:N - 1
    disp(i)

    % Runge-Kutta (2-nd order)
    k1 = h / 2 * (EvalPoly(Ht, xt_slave(i, :), sigma) + sync_coef .* (x(i, :) - xt_slave(i, :)));
    k2 = h * (EvalPoly(Ht, xt_slave(i, :) + k1, sigma) + sync_coef .* (x(i, :) - (xt_slave(i, :) + k1)));
    xt_slave(i + 1, :) = xt_slave(i, :) + k2;

    k1 = h / 2 * (EvalPoly(Ho, xo_slave(i, :), sigma) + sync_coef .* (x(i, :) - xo_slave(i, :)));
    k2 = h * (EvalPoly(Ho, xo_slave(i, :) + k1, sigma) + sync_coef .* (x(i, :) - (xo_slave(i, :) + k1)));
    xo_slave(i + 1, :) = xo_slave(i, :) + k2;
end

figure(1)
semilogy(t, vecnorm(x - xt_slave, 2, 2));
grid on;
title('Sync error (regular)')

figure(2)
semilogy(t, vecnorm(x - xo_slave, 2, 2));
grid on;
title('Sync error (orth)')