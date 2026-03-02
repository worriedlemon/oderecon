function H = wsindy(t, x, sigma, lambda)
% WSINDy - Weak Sparse Identification of Nonlinear Dynamics
% Inputs:
%   t      - [n x 1] time vector
%   x      - [n x d] state matrix (d = system dimension)
%   sigma  - [L x d] matrix of monomial powers for basis functions
%   M - number of test functions
%   lambda - sparsification parameter (threshold)
%   alpha - L2 regularization parameter
% Output:
%   H      - [L x d] matrix of identified coefficients
    
    h = t(2) - t(1);
    [n, vc] = size(x);
    M = size(sigma, 1);
    
    % Weak parameters
    width_sec = 0.5;                    % Width in seconds (fixed scale)
    width_pts = round(width_sec/h);     % Width in points
    step_pts = round(0.1*width_pts);    % Step is 10% of width
    centers = width_pts:step_pts:n-width_pts;
    K = length(centers);
    
    % Spline functions
    phi = @(s)(s>=0 & s<1).*(0.5*s.^2) + (s>=1 & s<2).*(0.75-(s-1.5).^2) + (s>=2 & s<=3).*(0.5*(3-s).^2);
    dphi = @(s)(s>=0 & s<1).*s + (s>=1 & s<2).*(3-2*s) + (s>=2 & s<=3).*(s-3);
    
    % Library of candidate functions
    Theta = zeros(n, M);
    for i = 1:M
        term = ones(n, 1);
        for j = 1:vc
            term = term .* (x(:,j).^sigma(i,j));
        end
        Theta(:,i) = term;
    end
    
    % Computation with normalization
    G = zeros(K,M); B = zeros(K,vc);
    
    for l = 1:K
        half_width_pts = ceil(width_pts/2);
        idx = (centers(l)- half_width_pts + 1):(centers(l) + half_width_pts);
        t_win = t(idx) - t(centers(l));              % Time relative to the center
        s_vals = 3*(t_win + width_sec/2)/width_sec;  % Scaling into [0,3]
        
        norm_coef = width_sec/3; % нормировка интеграла
        phi_vals = norm_coef*phi(s_vals); 
        dphi_vals = dphi(s_vals);
        
        % Интегралы с нормировкой
        for k = 1:vc
            B(l,k) = -sum(x(idx, k)' .* dphi_vals')*h; %Euler
            %B(l,k) = -trapz(x(idx, k)' .* dphi_vals')*dt; %Trapz
        end
        for j = 1:M
            G(l,j) = sum(Theta(idx,j) .* phi_vals)*h; %Euler
            %G(l,j) = trapz(Theta(idx,j) .* phi_vals)*dt; %Trapz
        end
    end
    
    % Original SINDy sparsification
    H = zeros(M,vc);
    for k = 1:vc
        Htemp = G \ B(:,k);
        for iter = 1:10
            small = abs(Htemp) < lambda;
            if (all(Htemp(small) == 0, "all"))
                break
            end

            Htemp(small) = 0;
            big = ~small;
            if sum(big) > 0
                Htemp(big) = G(:,big) \ B(:,k);
            end
        end
        H(:,k) = Htemp;
    end

end