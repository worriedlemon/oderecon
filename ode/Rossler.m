function dX = Rossler(t,X)
dX = X; 
a = 0.2;
b = 0.2;
c = 5.7;

x = X(1,:); y = X(2,:); z = X(3,:);

dX(1,:) = - y - z;
dX(2,:) = x + a*y;
dX(3,:) = b + z.*(x - c);
end