function varargout = prettyOrth(H,T,F,sigma,incf)
    % -- prettyorth(H, T, F, sigma)
    % -- prettyOrth(H, T, F, sigma, incf)
    % -- eqs = prettyorth(____)
    %     Prints out the equations for orthogonal
    %     polynomials, if the function called as procedure,
    %     otherwise saves equation strings in `eqs` cell.
    %
    %     H - coefficients cell
    %     T - degrees cell
    %     F - relation matrix, given by `orthpoly` function
    %     sigma - order ideal
    %     incf - logical value, indicating whether each
    %       equation should have 'f_n' prefix for n-th
    %       function
    
    [~, M] = size(H);
    [L, L] = size(F);
    
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
            str = char([]);
        end
        
        h = H{1,i};
        t = T{1,i};
        [~, vc] = size(t);
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
            
            [~, ind] = max(prod(sigma == t(j,:), 2));
            
            flag2 = sum(F(ind,:) ~= 0) > 1; 
            if flag2
                if flag1
                    str = [str, ' ('];
                    flag1 = 0;
                else
                    str = [str, '*('];
                end
            end
            
            bufstr = char([]);
            for k = 1:L % display entries
                if F(ind,k) ~= 0
                    if ~isempty(bufstr)
                        if (F(ind, k) > 0)
                            if F(ind,k) == 1
                                bufstr = [bufstr, ' +'];
                                flag1 = 1;
                            else
                                bufstr = [bufstr, ' + ', num2str(F(ind, k))];
                            end
                        else
                            if F(ind,k) == -1
                                bufstr = [bufstr, ' -'];
                                flag1 = 1;
                            else
                                bufstr = [bufstr, ' - ', num2str(-F(ind, k))];
                            end
                        end
                    elseif F(ind,k) < 0
                        if F(ind, k) == -1
                            bufstr = [bufstr, '-'];
                            flag1 = 1;
                        else
                            bufstr = [bufstr, '-', num2str(-F(ind, k))];
                        end
                    else
                        if flag2
                            bufstr = [bufstr, num2str(F(ind, k))];
                        else
                            bufstr = [bufstr, '*', num2str(F(ind, k))];
                        end
                    end
                    
                    for k2 = 1:vc
                        if sigma(k,k2) ~= 0
                            if flag1
                                bufstr = [bufstr, ' x', num2str(k2)]; %x1, x2, x3 ...
                                flag1 = 0;
                            else
                                bufstr = [bufstr, '*x', num2str(k2)]; %x1, x2, x3 ...
                            end
                            if sigma(k,k2) ~= 1
                                bufstr = [bufstr, '^', num2str(sigma(k,k2))];% x1^2, x2^3
                            end
                        end
                    end
                end
            end
            if isempty(bufstr) || strcmp(bufstr, '-')
                bufstr = [bufstr, '1'];
            end
            str = [str, bufstr];
            
            if flag2
                str = [str, ')'];
                flag2 = 0;
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
