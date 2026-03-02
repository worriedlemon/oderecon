function y = integrate_8ord(datarow, y0, htime)
%INTEGRATE_8ORD integrates datarow from initial state y0 over time htime
%   using 8th order Newton-Cotes method
%
%   INTEGRATE_8ORD(datarow, y0, h) - constant step size h
%   INTEGRATE_8ORD(datarow, y0, timespan) - variable step size
%
%   For the first few points where 8th order formula cannot be applied,
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
    % 0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 - 8
    c = [1070017, 4467094, -4604594, 5595358, -5033120, 3146338, -1291214, 312874, -33953]; 
    y(2) = y(1) + h*sum(c.*x(1:9)) / 3628800; % 0-1
    
    c = [-33953, 1375594, 3244786, -1752542, 1317280, -755042, 294286, -68906, 7297];
    y(3) = y(2) + h*sum(c.*x(1:9))  / 3628800; % 1-2
    
    c = [7297, -99626, 1638286, 2631838, -833120, 397858, -142094, 31594, -3233] ;
    y(4) = y(3) + h*sum(c.*x(1:9)) / 3628800; % 2-3

    c = [-3233, 36394, -216014, 1909858, 2224480, -425762, 126286, -25706, 2497] ;
    y(5) = y(4) + h*sum(c.*x(1:9))/ 3628800; % 3-4
    
    c = [2497, -25706, 126286, -425762, 2224480, 1909858, -216014, 36394, -3233] ;

    for j = 6:N-3
        xspan = x(j - 5:j + 3);
        y(j) = y(j-1) + h*sum(c.*xspan) / 3628800; % 4-5
    end

    c = [-3233, 31594, -142094, 397858, -833120, 2631838, 1638286, -99626, 7297] ;
    y(N-2) = y(N-3) + h*sum(c.*x(N-8:N)) / 3628800; % 5-6

    c = [7297, -68906, 294286, -755042, 1317280, -1752542, 3244786, 1375594, -33953];
    y(N-1) = y(N-2) + h*sum(c.*x(N-8:N))  / 3628800; % 6-7

    c = [-33953, 312874, -1291214, 3146338, -5033120, 5595358, -4604594, 4467094, 1070017];
    y(N) = y(N-1) + h*sum(c.*x(N-8:N))  / 3628800; % 7-8
end