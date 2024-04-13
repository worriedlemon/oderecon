function varargout = prettyABM(H,T,incf)
%PRETTYABM prints results of ABM - LSM reconstruction in a simple
%print-like form
%   PRETTYABM(H,T)
%   H is a cell array of size 1 x M, where M is dimension, containing lists
%   of coefficients by terms
%   T is a cell array of size 1 x M, where M is dimension, containing basis
%   terms obtained by delMinorTerms
%   incf is logical: if True, then prefix 'f_n = ' added to n_th equation
    [~, M] = size(H);
    
    if ~exist('incf', 'var')
        incf = 1;
    end
    
    if nargout > 0
        eqs = cell(1,M);
    end
    
    
    
    for i = 1:M %loop by number of functions
        if incf
            str = ['f_' , num2str(i), ' = ']; %string for entries
        else
            str = [];
        end
        
        h = H{1,i};
        t = T{1,i};
        [N, ~] = size(h); %number of terms
        for j = 1:N %loop by number of terms
            flag1 = 0; %flag, if 1 then do not draw * before monomial
            if j > 1 %formatting of the line
                if h(j) > 0
                    if ~isequal(num2str(h(j)),'1')
                        str = [str, ' + ', num2str(h(j))];
                    else
                        str = [str, ' +'];
                        flag1 = 1;
                    end
                else
                    if ~isequal(num2str(h(j)),'-1')
                        str = [str, ' - ', num2str(-h(j))];
                    else
                        str = [str, ' -'];
                        flag1 = 1;
                    end
                end
            else
                if isequal(num2str(h(j)),'-1')
                    str = [str, '-'];
                    flag1 = 1;
                else
                    if ~isequal(num2str(h(j)),'1')
                        str = [str, num2str(h(j))];
                    else
                        flag1 = 1;
                    end
                end
            end
            for k = 1:M %in each term, display entries
                if t(j,k) ~= 0
                    if flag1
                        str = [str, ' x', num2str(k)]; %x1, x2, x3 ...
                        flag1 = 0;
                    else
                        str = [str, '*x', num2str(k)]; %x1, x2, x3 ...
                    end
                    if t(j,k) ~= 1
                        str = [str, '^', num2str(t(j,k))];% x1^2, x2^3
                    end
                end
            end
        end
        
        if ~nargout
            disp(str);
        else
            eqs{1,i} = str;
        end
    end
    
    if nargout > 0
        varargout{1} = eqs;
    end
end
