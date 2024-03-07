% Used as a workaround for RNG to work in both Octave and Matlab
function rng_i(opt)
    if getenv('OCTAVE') == '1'
        switch opt
            case 'default'
                seed = 0;
            case 'shuffle'
                clk = clock();
                seed = clk(6);
            otherwise
                error('Argument could be "default" or "shuffle"');
        end
        rand('seed', seed);
    else
        rng opt;
    end
end
