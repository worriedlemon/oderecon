function dX = Duffing(t,X)
dX = X;
x = X(1,:);
y = X(2,:);

alpha = -1;
beta = 1;
delta = 0.02;
gamma = 3;
omega = 1;

dX(1,:) = y;
dX(2,:) = gamma * cos(omega * t) - (delta * y + alpha * x + beta * x ^ 3);
end