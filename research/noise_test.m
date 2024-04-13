rng_i default;
warning off;

system = @Lorenz;
eqs = 1;
calc = 0;
drawplt = 0;

Tmax = 100; % Time end
h = 0.01; % Step
start_point = [4 -2 0]; % Initial point

timespan = 0:h:Tmax;
[~, x] = ode45(system, timespan, start_point);
y = transpose(system(0, x'));
#y = diff4(x, timespan);

vc = length(start_point);
deg = 2;

noise_amp = input(['Input desired noise magnitude:', newline, '> ']);

N = 12;
ri = randi(size(y, 1), 1, N);
rx = x(ri,:) + noise_amp * randn(N, vc);
ry = y(ri,:) + noise_amp * randn(N, 1);
clear ri noise_amp;

polynomial_options = {'x', 'bernstein', 'orth'};
fnum = listdlg("PromptString", 'Choose desired polynomial basis', "ListString", polynomial_options);
assert(~isempty(fnum), 'Program execution stopped due to no given answer');

opt = {polynomial_options{1, fnum}};
disp(["Using \"", opt{1,1}, "\" polynomials"]);

eta = 1e-5;
switch opt{1,1}
    case 'x'
        tx = rx;
        [~, sigma] = ApproxBM(x, eta, deglexord(deg, vc));
    case 'bernstein'
        c01 = repmat([0; 1], 1, vc);
        [tx, c] = affine_transform(rx, c01);
        sigma = berndeg(deg, vc);
    case 'orth'
        tx = rx;
        [~, sigma] = ApproxBM(x, eta, deglexord(deg, vc));
        
        orthinterv = [0, 1];
        [F, nrms] = orthpoly(deg, vc, orthinterv(1), orthinterv(2), 1e-9, 0);
        opt = {opt{1,1}, F, sigma};
        orthogonality_test(opt{1,2}, deg, vc, orthinterv(1), orthinterv(2), 'brief', 1e-6, nrms);
        clear orthinterv F;
end
clear eta polynomial_options;

try
    close(fnum);
end

H = cell(1, vc);
T = cell(1, vc);
tol = 1e-5;

disp('Reconstruction started...')
for i = 1:vc
    [hi, tau] = delMinorTerms(tx, ry(:,i), sigma, tol, opt);
    norm(ry(:,i) - EvalPoly(hi, tx, tau, opt))
    H{1,i} = hi;
    T{1,i} = tau;
end
disp("OK.");

if calc
    disp("\nFinding values based on reconstructed equation...")
    switch opt{1,1}
        case 'x'
            [~, x1] = ode45(@(t,x)oderecon(H, T, t, x, opt), timespan, start_point);
        case 'orth'
            [~, x1] = ode45(@(t,x)oderecon(H, T, t, x, opt), timespan, start_point);
        case 'bernstein'
            [~, x1] = ode45(@(t,x)oderecon(H, T, t, affine_transform(x', c01, c)', opt), timespan, start_point);
    end
end

if eqs
    disp("\nEquations:")
    switch opt{1,1}
        case 'x'
            prettyABM(H, T);
        case 'orth'
            prettyOrth(H, T, opt{1,2}, sigma);
        case 'bernstein'
            prettyBernstein(H, T);
    end
end

if drawplt && calc
    figure(fnum);
    plot3(x(:,1), x(:,2), x(:,3), 'b', start_point(1), start_point(2), start_point(3), 'r*');
    hold on;
    scatter3(rx(:,1), rx(:,2), rx(:,3), 'k');
    xlabel('\itx'); ylabel('\ity'); zlabel('\itz');
    plot3(x1(:,1), x1(:,2), x1(:,3), 'y');
    legend('Initial dynamic system', 'Start point', 'Random points', 'Reconstructed system');
    title(sprintf('System reconstruction (%s)', opt{1,1}));
 end