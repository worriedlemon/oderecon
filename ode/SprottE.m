function dX = SprottE(t,X)
  %initial can be [0.5 0.2 1]
  dX = X; 
  a = 1.01;

  x = X(1,:); y = X(2,:); z = X(3,:);

  dX(1,:) = a*y.*z;
  dX(2,:) = x.^2 - y;
  dX(3,:) = 1 - 4*x;
end
