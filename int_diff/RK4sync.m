function Ys = RK4sync(Y,H,T,Ts,synch_coef)
%Solve synchronized system by RK4 method, with Lagrangian interpolation
%between points of signal in times 1/2
%input:
% Y - data of master system
% H,T - data of the reconstructed system via delMinorTerms method
% Ts - discretization of Y
%output:
% Ys - solution of ODE


%k_sync = [synch_coef,synch_coef,synch_coef];
k_sync = [synch_coef,0,0];

[N,~] = size(Y);

Ys = zeros(N,3);
Ys(1,:) = Y(1,:);
U = Ys(1,:);

H1 = H{1,1};
H2 = H{1,2};
H3 = H{1,3};

T1 = T{1,1};
T2 = T{1,2};
T3 = T{1,3};

ctr = 2;
for i = 2:N
    h = Ts(i) - Ts(i-1); %time step
    if (i <= 2) || (i >= N - 1)
        Y12 = 0.5*(Y(i,:) + Y(i-1,:)); %first-order hold of signal
    else
        %Lagrange polynomial
        x = Ts(i-2:i+1);
        y = Y(i-2:i+1,:);
        sum = 0;
        a = Ts(i-1) + 0.5*h;
        for ic = 1:length(x)
            u = 1;
            l = [1,1,1];
            for j = 1:length(x)
                if j ~= ic
                    u = u * (a - x(j));
                    l = l.* (x(ic) - x(j));
                end
            end
            sum= sum + u ./ l.* y(ic,:);
        end       
        Y12 = sum;
    end

%     %proportional control
    k1 = [EvalPoly(H1,U,T1),EvalPoly(H2,U,T2),EvalPoly(H3,U,T3)] - k_sync.*(U - Y(i-1,:));

    U1 = U + h/2*k1;
    k2 = [EvalPoly(H1,U1,T1),EvalPoly(H2,U1,T2),EvalPoly(H3,U1,T3)] - k_sync.*(U1 - Y12);

    U2 = U + h/2*k2;
    k3 = [EvalPoly(H1,U2,T1),EvalPoly(H2,U2,T2),EvalPoly(H3,U2,T3)] - k_sync.*(U2 - Y12);

    U3 = U + h*k3;
    k4 = [EvalPoly(H1,U3,T1) ,EvalPoly(H2,U3,T2) ,EvalPoly(H3,U3,T3)]  - k_sync.*(U3 - Y(i,:));

    %Pecora-Caroll control
    % nv = 1;
    % 
    % U(nv) = Y(i-1,nv);
    % 
    % k1 = [EvalPoly(H1,U,T1),EvalPoly(H2,U,T2),EvalPoly(H3,U,T3)];
    % 
    % U1 = U + h/2*k1;
    % U1(nv) = Y12(nv);
    % 
    % k2 = [EvalPoly(H1,U1,T1),EvalPoly(H2,U1,T2),EvalPoly(H3,U1,T3)];
    % 
    % U2 = U + h/2*k2;
    % U2(nv) = Y12(nv);
    % 
    % k3 = [EvalPoly(H1,U2,T1),EvalPoly(H2,U2,T2),EvalPoly(H3,U2,T3)];
    % 
    % U3 = U + h*k3;
    % U3(nv) = Y(i,nv);
    % 
    % k4 = [EvalPoly(H1,U3,T1) ,EvalPoly(H2,U3,T2) ,EvalPoly(H3,U3,T3)];
    
    U = U + h/6*(k1 + 2*k2 + 2*k3 + k4);
    Ys(ctr,:) = U;
    ctr = ctr + 1;
end