function [sc_x, from_a, from_b] = affine_transform(varargin)
  % Syntax:
  % sc_x = affine_transform(x, to_a, to_b)
  % sc_x = affine_transform(x, to_a, to_b, from_a, from_b)
  % [sc_x, from_a, from_b] = affine_transform(____)
  
  if ~(nargin == 3 || nargin == 5)
    error(sprintf("Function receives 3 or 5 arguments (%u given)", nargin))
  end
  
  x = varargin{1,1};
  assert(size(x, 1) == 1 || size(x, 2) == 1, "First argument should be a 1d-vector");
  to_a = varargin{1,2};
  to_b = varargin{1,3};
  
  if nargin == 5
    from_a = varargin{1,4};
    from_b = varargin{1,5};
  else
    from_a = min(x);
    from_b = max(x);
  end
  
  sc_x = to_a + (to_b - to_a) / (from_b - from_a) * (x - from_a);
end