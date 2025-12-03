function dX = Cang(t,X)
dX = X; 
a = 50;
b = -16;
c = 10;
d = 0.2;
e = 10;
f = 16;
g = 0.5;

x = X(1,:); y = X(2,:); z = X(3,:); w = X(4,:);

dX(1,:) = -a * x - e * w + y .* z;
dX(2,:) = b * y + x .* z;
dX(3,:) = c * z + f * w - x .* y;
dX(4,:) = d * w - g * z;
end