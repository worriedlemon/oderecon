function setup(varargin)
    disp("Current work directory is:");
    disp(pwd);
    disp("");
    
    paths = {"examples", "misc", "ode", "int_diff"};
    
    for i = 1:size(paths, 2)
        incl = strcat(pwd, filesep, paths{1, i});
        addpath (incl);
        disp(strcat("Path added :", incl));
    end
    
    if nargin == 1 && (toupper (varargin{1,1})) == "OCTAVE"
        disp("\nAssuming you are using Octave");
        setenv("OCTAVE", "1");
        disp("New environment variable added: OCTAVE=1");
    endif
    
    disp("");
end
