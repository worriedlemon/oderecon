function dX = Thomas(t,X)
  %initial can be [0.5 0.2 1]
  dX = X; 
  b = 0.185;

  x = X(1,:); y = X(2,:); z = X(3,:);

  dX(1,:) = sin(y) - b*x;
  dX(2,:) = sin(z) - b*y;
  dX(3,:) = sin(x) - b*z;
end