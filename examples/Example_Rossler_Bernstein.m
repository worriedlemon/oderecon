% Initializing Random Number Generator (octave support)
rand_init;
close all;

% Rossler system Simulation
Tmax = 100; % Time end
h = 0.01; % Step
start_point = [4, -2, 0]; % Initial point
disp("Start point is "), disp(start_point)

[t,y] = ode45(@Rossler, [0:h:Tmax], start_point); % Getting attractor points array via Runge-Kutta ODE solver
w = transpose(Rossler(0, y')); % Getting differentials for each point
[y, ymin, ymax] = normalizeTo(y, 0, 1); 

% x-z plane projection
figure(1);
plot3(y(:,1),y(:,2),y(:,3));
xlabel('\itx');
ylabel('\ity');
zlabel('\itz');

N = 15; % Data points
[count, M] = size(y); % Points count and dimensions
W = zeros(N,M);
Y = zeros(N,M);

for i = 1:N % Taking random points from attractor
    id = ceil(rand * count);  % Index of data point
    W(i,:) = w(id,:); % X
    Y(i,:) = y(id,:); % Y
end


% Plot sample points
figure(1); hold on
scatter3(Y(:,1),Y(:,2),Y(:,3),23,'MarkerEdgeColor','g','MarkerFaceColor','y','LineWidth',1.5);

% Reconstruct order ideal
eps = 1e-5;
[~, O] = ApproxBM(Y, eps, deglexord(2,3)) % Use approximate Buchberger-Moller algorithm

eta = 1e-7;
H = cell(1,3);
T = cell(1,3);

error("Not implemented yet")

%reconstruct each equation
for i = 1:3
    V = W(:,i);
    #[hi,tau] = delMinorTerms(Y,V,O,eta); %get equation and basis
    V0 = EvalPoly(hi,Y,tau); 
    norm(V - V0) %check if norm is appropriate
    
    H{1,i} = hi;
    T{1,i} = tau;
end

%simulate results
[~,y] = ode45(@(t,x)oderecon(H,T,t,x),[0:h:Tmax], start_point); %solve ODE
figure(1);
plot3(y(:,1),y(:,2),y(:,3),'-');

%display equations
prettyABM(H,T)