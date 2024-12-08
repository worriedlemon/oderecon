function setup(varargin)
    % -- setup args...
    %     Used to set oderecon project up before using scripts.
    %     
    %     Arguments can be:
    %     path - adds all directories to PATH environment variable
    %     octave - sets environment variable OCTAVE=1 (is used as
    %       convenience method for functions to run in GNU Octave)
    %     format - sets format to 'short g'
    
    if (nargin == 0)
      disp('Possible arguments: path, octave, format');
      return
    end
    
    argument = lower(varargin{1,1});
    
    if strcmp(argument, 'path')
        disp('Current working directory is:');
        disp(pwd);
        disp('');
        
        dirl = dir;
        paths = {''};
        for i = 1:length(dir)
            ignore = length(dirl(i).name) >= 4 && strcmp(dirl(i).name(1:4), 'IGN_');
            github = length(dirl(i).name) >= 7 && strcmp(dirl(i).name(1:7), 'GITHUB_');
            if isfolder(dirl(i).name) && dirl(i).name(1) ~= '.' && ~ignore && ~github
                paths = {paths{1:length(paths)}, dirl(i).name};
            end
        end
        
        for i = 1:size(paths, 2)
            incl = strcat(pwd, filesep, paths{1, i});
            addpath(incl);
            disp(strcat('Path added :', incl));
        end
    elseif strcmp(argument, 'octave')
        setenv('OCTAVE', '1');
        disp('Assuming you are using Octave: new environment variable added (OCTAVE=1)');
    elseif strcmp(argument, 'format')
        format short g;
        disp('Set "format short g"');
    else
        warning('Invalid argument "%s" skipped', argument);
    end
    
    disp('');
    
    if nargin > 1
        setup(varargin{1,2:nargin});
    end
end
