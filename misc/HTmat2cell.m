function [Hcell, Tcell] = HTmat2cell(H, T)
    [mc, eqc] = size(T);
    Tcell = mat2cell(repmat(T, 1, eqc), mc, repmat(eqc, 1, eqc));
    Hcell = mat2cell(H, mc, ones(1, eqc));
end

