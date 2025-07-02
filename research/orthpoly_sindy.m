rng_i default;
warning off;
close all

% System coefficients
Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.3 0 0 0 0 0 0 0; 0.3 0 0 -5.7 0 0 1 0 0 0]';

% Used system
%system = @Rossler
system = @Lorenz

Href = eval(['H', lower(func2str(system))]); % Used coefficients
Ttrans = 100; % Transient time
Tmax = 5; % Time of experimental series for Lorenz
%Tmax = 10; % Time of experimental series for Rossler
h = 1e-2; % Step
start_point = [4 -2 0]; % Initial point

% %first, go transient
[~, x] = ode78(system, 0:h:Ttrans, start_point);
start_point = x(end,:);
%then, go experiment
[t, x] = ode78(system, 0:h:Tmax, start_point);

y = transpose(system(0, x'));

[N, vc] = size(x);

eqc = size(y, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

% initial solution
F = orthpoly_t(sigma, t, x);
Ho = zeros(mc, eqc);
E = EvalPoly(F', x, sigma);
for i = 1:eqc
    Ho(:, i) = trapz(x(:, i), E);
    % for j = 1:mc
    %     Ho(:, i) = intdiff4(x(:, i), E);
    % end
end

% using SINDy approach (no derivatives), improving the solution
H = Ho;
lambda = 1;
k = 0;
while (k < 10)
    smallinds = (abs(H) < lambda);
    % all small values are already zeros, stop
    if (all(H(smallinds) == 0, "all"))
        break
    end

    H(smallinds) = 0;
    Ht = zeros(mc, vc);
    for ind = 1:vc
        biginds = ~smallinds(:, ind);
        sigma_temp = sigma(biginds, :);
        F_temp = orthpoly_t(sigma_temp, t, x);
        E = EvalPoly(F_temp', x, sigma_temp);
        H(biginds, ind) = trapz(x(:, ind), E);
        Ht(biginds, ind) = F_temp' * H(biginds, ind);
    end
    k = k + 1;
end

disp("Result: ");
disp(Ht);