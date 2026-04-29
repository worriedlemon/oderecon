function x = fcompose(fun, n, x)
    % FCOMPOSE - compose function fun n times with x argument

    assert(n >= 0, 'n should be greater or equal 1')
    
    for i = 1:n
        x = fun(x);
    end
end

