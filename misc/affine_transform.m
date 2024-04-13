function [sc_x, c_from] = affine_transform(x, c_to, c_from)
  % -- sc_x = affine_transform(x, c_to)
  % -- sc_x = affine_transform(x, c_to, c_from)
  % -- [sc_x, c_from] = affine_transform(____)
  %
  %     x is input matrix of size N x M
  %     c_to, c_from is domain matrices of size 2 x M (first row contains lower bound, second row contains higher bound)
  %     sc_x is scaled matrix of size N x M
    if ~exist('c_from', 'var')
        c_from = [min(x); max(x)];
    end
  
    sc_x = c_to(1,:) + (c_to(2,:) - c_to(1,:)) ./ (c_from(2,:) - c_from(1,:)) .* (x - c_from(1,:));
end
