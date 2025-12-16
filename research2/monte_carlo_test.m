rng_i shuffle % maybe default is better
warning off
close all

% Used system
sys = @Rossler

sysname = func2str(sys);

[Href, deg, vc, mc] = load_href(sysname); % Used coefficients

T_range = [100; 150];    % relatively large T
h_range = [1e-3; 5e-2];  % relatively small h
x_range = [-5; 5];
c01 = [0; 1];

sigma = deglexord(deg, vc);

addpath('research2/sindy_interface/')

% Methods:
% 1. SINDY
% 2. Orthpoly-SINDY
% 3. Ohrtpoly-Int-SINDY
methods = get_methods;
mtdn = size(methods, 2);

% Lambda parameters
calculate_lambdas = 0;
default_lambdas = 1;
if (calculate_lambdas)
    disp('Calculating lambdas:')
    lambdas = get_lambdas(methods, sys, sigma, Href);
elseif (default_lambdas)
    default_lambda = 1e-2; % Here is the 
    disp('Using default lambdas:')
    lambdas = default_lambda * ones(1, mtdn) %#ok
else
    lambdas_rossler = [5.29831690628370e-5, 0.000573615251044868, 0.000174332882219999];
    lambdas_lorenz = [0.0137382379588326, 0.0137382379588326, 0.100000000000000];
    
    disp('Using precalculated lambdas:')
    lambdas = eval(['lambdas_', lower(sysname)]) %#ok
end

% Pairwise compared methods
methods_pairs = nchoosek(1:mtdn, 2);
methods_pairs_size = size(methods_pairs, 1);

tests = 5;
pair_resultsXi = -ones(mtdn);
pair_resultsR = pair_resultsXi;
for p = 1:methods_pairs_size
    [first_method, second_method] = methods{1, methods_pairs(p, :)};
    disp(['Test method pair ', num2str(p), ': ', method2str(first_method), ' vs ', method2str(second_method)])
    
    winsXi = 0;
    winsR = 0;
    for k = 1:tests
        disp(['Test #', num2str(k)])
        
        % Take random
        Tmax = affine_transform(rand(1), T_range, c01);
        start_point = affine_transform(rand(1, vc), x_range, c01);
        h = affine_transform(rand(1), h_range, c01);
        [t, x_ref] = ode78(sys, 0:h:Tmax, start_point);
    
        [Xi1, Rscript1] = get_norms(Href, t, x_ref, sigma, lambdas(methods_pairs(p, 1)), first_method);
        [Xi2, Rscript2] = get_norms(Href, t, x_ref, sigma, lambdas(methods_pairs(p, 2)), second_method);
    
        winsXi = winsXi + (Xi1 < Xi2);
        winsR = winsR + (Rscript1 < Rscript2);
    end

    pair_resultsXi(methods_pairs(p, 1), methods_pairs(p, 2)) = winsXi;
    pair_resultsXi(methods_pairs(p, 2), methods_pairs(p, 1)) = tests - winsXi;
    pair_resultsR(methods_pairs(p, 1), methods_pairs(p, 2)) = winsR;
    pair_resultsR(methods_pairs(p, 2), methods_pairs(p, 1)) = tests - winsR;
end

disp('Results matrixes:');
disp('Xi: '); disp(cat(2, pair_resultsXi, sum(pair_resultsXi, 2) + 1))
disp('Rho: '); disp(cat(2, pair_resultsR, sum(pair_resultsR, 2) + 1))

rmpath('research2/sindy_interface/')