function setup(varargin)
    if (nargin == 0)
      disp('Possible arguments: path, octave, format');
      return
    end
    
    argument = lower(varargin{1,1});
    
    if strcmp(argument, 'path')
        disp('Current working directory is:');
        disp(pwd);
        disp('');
        
        paths = {'', 'examples', 'misc', 'ode', 'int_diff', 'orthpoly'};
        
        for i = 1:size(paths, 2)
            incl = strcat(pwd, filesep, paths{1, i});
            addpath (incl);
            disp(strcat('Path added :', incl));
        end
    elseif strcmp(argument, 'octave')
        setenv('OCTAVE', '1');
        disp('Assuming you are using Octave: new environment variable added (OCTAVE=1)');
    elseif strcmp(argument, 'format')
        format short g;
        disp('Set "format short g"');
    else
        warning(sprintf('Invalid argument "%s" skipped', argument));
    end
    
    disp('');
    
    if nargin > 1
        setup(varargin{1,2:nargin});
    end
end
