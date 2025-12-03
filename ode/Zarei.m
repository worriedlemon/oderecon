function dX = Zarei(t,X)
dX = X; 
a = 10;
b = 60;
c = 20;
d = 15;
e = 40;
f = 1;
g = 50;
h = 10;

x = X(1, :);
y = X(2, :);
z = X(3, :);
u = X(4, :);
v = X(5, :);

dX(1, :) = -a * x + y .* z;
dX(2, :) = -b * y + f * v;
dX(3, :) = -c * z + g * u + x .* y;
dX(4, :) = d * u - h * z;
dX(5, :) = e * v - y .* x .* x;
end