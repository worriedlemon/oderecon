function varargout = prettyBernstein(H, T, varargin)
    [~, M] = size(H);
    str = [];
    incf = 1;
    delim = ';';
    nargin = nargin - 2;
    
    if nargin >= 1
        incf = varargin{1,1};
    end
    
    if nargin >= 2
        delim = varargin{1,2};
    end
    
    for i = 1:M %loop by number of functions
        if incf
            str = ['f_' , num2str(i), ' = ']; %string for entries
        end
        
        h = H{1,i};
        t = T{1,i};
        [N, ~] = size(h); %number of terms
        for j = 1:N %loop by number of terms
            flag1 = 0; %flag, if 1 then do not draw * before monomial
            if j > 1 %formatting of the line
                if h(j) > 0
                    if h(j) ~= 1
                        str = [str, ' + ', num2str(h(j))];
                    else
                        str = [str, ' +'];
                        flag1 = 1;
                    end
                else
                    if h(j) ~= -1
                        str = [str, ' - ', num2str(-h(j))];
                    else
                        str = [str, ' -'];
                        flag1 = 1;
                    end
                end
            else
                if h(j) == -1
                    str = [str, '-'];
                    flag1 = 1;
                else
                    if h(j) == 1
                        str = [str, num2str(h(j))];
                    else
                        flag1 = 1;
                    end
                end
            end
            
            % Calculating binomial coefficient
            binom = 1;
            for k = 1:M
                binom = binom * nchoosek(t(j, k), t(j, k + M));
            end
            if binom ~= 1
                if flag1
                    str = [str, num2str(binom), '*'];
                    flag1 = 0;
                else
                    str = [str, '*', num2str(binom), '*'];
                end
            elseif ~flag1
                str = [str, '*'];
            end
            
            % Creating a monomial from each term
            buf = '';
            for k = 1:M
                buf2 = '';
                nk = t(j, k); ik = t(j, k + M);
                if ik ~= 0
                    buf2 = [buf2, 't', num2str(k)];
                    if ik ~= 1
                        buf2 = [buf2, '^', num2str(ik)];
                    end
                end
                if nk ~= ik
                    if (ik ~= 0)
                        buf2 = [buf2, '*'];
                    end
                    buf2 = [buf2, '(1 - t', num2str(k), ')'];
                    if nk - ik ~= 1
                        buf2 = [buf2, '^', num2str(nk - ik)];
                    end
                end
                if ~isempty(buf) && ~isempty(buf2)
                   buf = [buf, '*'];
                end
                buf = [buf, buf2];
            end
            str = [str, buf];
        end
        
        if ~nargout
            disp(str);
            str = [];
        elseif i ~= M
            str = [str, delim];
        end
    end
    
    if nargout > 0
        varargout{1} = str;
    end
end