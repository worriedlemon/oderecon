rng_i default;
warning off;

system = @Lorenz;

Tmax = 50; % Time end
h = 0.01; % Step
start_point = [4 -2 0]; % Initial point

[t, x] = ode45(system, 0:h:Tmax, start_point);
y = diff4(x, t);

noise_amp = 0.001;

[N, vc] = size(x);

rx = x + noise_amp * randn(N, vc);
%ry = diff4(rx, t);
ry = y;

eqc = size(ry, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);
F = orthpoly_t(sigma, t, x);

H = zeros(mc, eqc);
tol = 1e-7;

for i = 1:eqc
    [H(:, i), ~] = delMinorTerms_dx(rx(:, i), rx, ry(:, i), F, sigma, tol, 0);
end

nrm = vecnorm(ry - EvalPolyOrth(H, rx, sigma, F, sigma))

H = F' * H