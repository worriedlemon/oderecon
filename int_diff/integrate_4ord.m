function y = integrate_4ord(datarow, y0, htime)
%INTEGRATE_4ORD integrates datarow from initial state y0 over time htime
%   using 4th order Newton-Cotes method (Boole's rule)
%
%   INTEGRATE_4ORD(datarow, y0, h) - constant step size h
%   INTEGRATE_4ORD(datarow, y0, timespan) - variable step size
%
%   For the first few points where 4th order formula cannot be applied,
%   method falls back to lower order formulas

     % Handle constant step size input
    if length(htime) == 1
        h = htime;
    else
        h = mean(diff(htime));
    end
    
    x = datarow(:)';  % Ensure row vector
    N = length(x);
    y = zeros(1, N);
    y(1) = y0;
    
    % For first few points, use progressively higher order methods
    % 0 - 1 - 2 - 3 - 4
    c = [251, 646, -264, 106, -19] / 720; 
    y(2) = y(1) + h*sum(c.*x(1:5)); % 0-1
    
    c = [-19, 346, 456, -74, 11] / 720;
    y(3) = y(2) + h*sum(c.*x(1:5)); % 1-2
    
    c = [11, -74, 456, 346, -19] / 720;
    for j = 4:N-1
        xspan = x(j - 3:j + 1);
        y(j) = y(j-1) + h*sum(c.*xspan); %2-3
    end

    c = [-19, 106, -264, 646, 251] / 720; 
    y(N) = y(N-1) + h*sum(c.*x(N-4:N)); % 3-4
end