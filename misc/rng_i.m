function rng_i(opt)
    % -- rng_i default
    % -- rng_i shuffle
    % -- rng_i(seed)
    %
    %     Argument is 'default' for setting RNG seed to 0
    %     Argument is 'shuffle' for setting RNG sead to system time
    %     Argument is integer value 'seed', RNG sets seed to specified value
    % 
    %     Note: used as a workaround for RNG to work in both Octave and Matlab
    str = isstring(opt) || ischar(opt);
    
    if getenv('OCTAVE') == '1'
        if ~str
            seed = opt;
        else
            switch opt
                case 'default'
                    seed = 0;
                case 'shuffle'
                    clk = clock();
                    seed = clk(6);
                otherwise
                    error('Argument could be "default", "shuffle" or seed value');
            end
        end
        rand('seed', seed);
        randn('seed', seed);
    else
        if ~str
            rng(opt);
        else
            rng opt;
        end
    end
end
