% Used as a workaround for RNG to work in both Octave and Matlab
function rand_init()
    if getenv("OCTAVE") == "1"
      disp('Assuming you using Octave');
      rand("seed", (clock())(6));
    elseif
      rng default;
    endif
end
