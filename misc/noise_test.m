rng_i default;
warning off;

system = @Lorenz;
Tmax = 100; % Time end
h = 1e-3; % Step
start_point = [4 -2 0]; % Initial point

timespan = 0:h:Tmax;
[~, x] = ode45(system, timespan, start_point);
y = transpose(system(0, x'));

fc = 3;
deg = 2;
c01 = repmat([0; 1], 1, fc);

noise_amp = input(['Input desired noise magnitude:', newline, '> ']);

N = 25;
ri = randi(size(y, 1), 1, N);
rx = x(ri,:) + noise_amp * randn(N, fc);
ry = y(ri,:) + noise_amp * randn(N, 1);
clear ri;

polynomial_options = {'x', 'bernstein'};
ans = listdlg("PromptString", 'Choose desired polynomial basis', "ListString", polynomial_options);
opt = polynomial_options{1,ans};
disp(["Using \"", opt, "\" polynomials"]);

switch opt
    case 'x'
        tx = rx;
        [~, sigma] = ApproxBM(x, 1e-5, deglexord(deg, fc));
    case 'bernstein'
        [tx, c] = affine_transform(rx, c01);
        sigma = berndeg(deg, fc);
end

clear polynomial_options;

try
    close(ans);
end
figure(ans);
plot3(x(:,1), x(:,2), x(:,3), 'b', start_point(1), start_point(2), start_point(3), 'r*');
hold on;
scatter3(rx(:,1), rx(:,2), rx(:,3), 'k');
xlabel('\itx'); ylabel('\ity'); zlabel('\itz');

H = cell(1, fc);
T = cell(1, fc);
eta = 1e-3;

disp('Reconstruction started...')
for i = 1:fc
    [hi, tau] = delMinorTerms(tx, ry(:,i), sigma, eta, opt);
    norm(ry(:,i) - EvalPoly(hi, tx, tau, opt))
    H{1,i} = hi;
    T{1,i} = tau;
end

disp("OK.\n\nFinding values based on reconstructed equation...")
switch opt
    case 'x'
        [~, x1] = ode45(@(t,x)oderecon(H, T, t, x, opt), 0:h:Tmax, start_point);
    case 'bernstein'
        [~, x1] = ode45(@(t,x)oderecon(H, T, t, affine_transform(x', c01, c)', opt), 0:h:Tmax, start_point);
end

disp("Done.\n\nAverage system reconstruction errors:")
err = mean(x1 - x)

plot3(x1(:,1), x1(:,2), x1(:,3), 'y');
legend('Initial dynamic system', 'Start point', 'Random points', 'Reconstructed system');
title(sprintf('System reconstruction (%s)\nAverage errors: %s', opt, disp(err)));