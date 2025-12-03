function [Pref, vc] = load_pref(name)
    Psprotte = [0.76, 0.22, 3.00];
    Prossler = [7.72, -10.67, 0.31];
    Plorenz = [4.08, -2.21, 30.2];
    Pcang = [3.84, 6.81, 24.06, -3.46];
    Pzarei = [138.98, 26.92, -107.96, -51.46, 12463.08];
    Pvdpl = [-2, 0];

    pname = ['P', lower(name)];
    Pref = eval(pname);
    vc = length(Pref);
end

