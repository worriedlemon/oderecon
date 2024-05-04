rng_i default;
close all;
warning off;

% Dynamic system Simulation
start_point = [4 -2 0]; % Initial point
h = 0.01; % Step
Tspan = [10, 100:100:1000] * h * 10;

deg = 2; % Degree of reconstructed function
vc = 3; % Variables count
eqc = 3; % Equations count
mc = nchoosek(deg + vc, vc); % Monomials count

sigma = deglexord(deg, vc);

Ns = zeros(length(Tspan), 1);
Es = [Ns, Ns];

for i = 1:length(Tspan) % Time max
    disp(i)
    Tmax = Tspan(i);
    [t, x] = ode45(@Rossler, 0:h:Tmax, start_point);
    y = transpose(Rossler(0, x'));
    Ns(i) = length(t);

    tic; % Timer start
    [F, norms] = orthpoly_t(sigma, t, x, 0);
    coefs = zeros(mc, eqc);
    for eq = 1:eqc
        for j = 1:mc
            coefs(j, eq) = trapz(x(:, eq), EvalPoly(transpose(F(j, :) / norms(j)), x, sigma));
        end
    end

    %coefs = F' * coefs
    Es(i, 1) = toc;
    
    tic;
    E = EvalPoly(eye(mc), x, sigma);
    coefs = (E'*E)\E'*y;
    Es(i, 2) = toc; % timer end
end

semilogy(Ns, Es(:, 1), 'b', Ns, Es(:, 2), 'r');
xlabel('Points count'); ylabel('Elapsed time, \t{s}');
title('Speed comparison between LSM and Orth');
legend('Orth', 'LSM');