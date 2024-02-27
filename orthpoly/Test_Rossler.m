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

figure(1, "Position", [200 200 960 720]);
plot3(x(:,1), x(:,2), x(:,3), "b", start_point(1), start_point(2), start_point(3), "r*");
hold on;
scatter3(rx(:,1), rx(:,2), rx(:,3), "k");
xlabel("\itx"); ylabel("\ity"); zlabel("\itz");
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
for i = 1:fc
  tx1(:,i) = affine_transform(x(:,i), 0, 1, c(1,i), c(2,i));
end
y1 = bernbase(tx1', sigma) * coefs;

disp("Reconstructed system error:");
err = norm(y - y1)

function dx = sda(x, sigma, coefs, fc, c)
  for i = 1:fc
    tx1(:,i) = affine_transform((x')(:,i), 0, 1, c(1,i), c(2,i));
  end
  dx = transpose(bernbase(tx1', sigma) * coefs);
end

[~, x1] = ode45(@(t, x)sda(x, sigma, coefs, fc, c), 0:h:Tmax, start_point);

plot3(x1(:,1), x1(:,2), x1(:,3), "-y", "DisplayName", "Reconstructed function");
title(sprintf("Rossler system reconstruction\nError: %g", err));
