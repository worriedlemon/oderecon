rand_init;
close all;

Tmax = 50;
h = 0.01;
start_point = [4 -2 0];
[t, x] = ode45(@Rossler, 0:h:Tmax, start_point);
y = transpose(Rossler(0, x'));

N = 25;
ri = randi(size(y, 1), 1, N);
rx = x(ri,:);
ry = y(ri,:);

figure(1);
plot3(x(:,1), x(:,2), x(:,3), "b", start_point(1), start_point(2), start_point(3), "r*");
hold on;
scatter3(rx(:,1), rx(:,2), rx(:,3), "k");
title("Function reconstruction");
legend("Rossler attractor", "Start point", "Random points");

fc = size(y, 2);
deg = 2;
sigma = berndeg(deg, fc);

tx = rx;
c = zeros(2, fc);
for i = 1:fc
  [tx(:,i), c(1,i), c(2,i)] = affine_transform(rx(:,i), 0, 1);
end

B = bernbase(tx', sigma);
coefs = (B'*B)\B'*ry

tx1 = x;
y1 = y;
for i = 1:fc
  tx1(:,i) = affine_transform(x(:,i), 0, 1, c(1,i), c(2,i));
  y1(:,i) = bernbase(tx1', sigma) * coefs(:,i);
end

norm(y - y1)