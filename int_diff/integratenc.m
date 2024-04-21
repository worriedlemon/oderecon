function yint = integratenc(t,y,n,ord)
%INTEGRATENC integrates with general Newton-Cotes formula
%   yint = INTEGRATENC(t,y,n,ord)
%   t - time span
%   y - y span
%   n - number of data points in a formula
%   ord - order of approximating polynomial

N = length(t); %number of data points
p = ord + 1; %length of formula (order of parabola + 1)

if p > N %if formula is longer than data span, correct it
    p = N;
end

yint = zeros(N,1);

n2 =  ceil((n + 1)/2);

%take first n/2 points

for j = 2:n2 
    tspan = t(1:n);
    yspan = y(1:n);
    istart = j - 1;
    iend = j;
    yint(j) = yint(j-1)+localnc(tspan,yspan,istart,iend,n,p);
end

%then, use sliding formula of length p

for j = n2 + 1: N - n2 %take p points
    tspan = t(j - n2 + 1:j - n2 + n);
    yspan = y(j - n2 + 1:j - n2 + n);
    istart = n2 - 1;
    iend = n2; %center
    yint(j) = yint(j-1)+localnc(tspan,yspan,istart,iend,n,p);
end

%then, take last p/2 points

for j = N - n2 + 1: N
    tspan = t(N - n + 1:N);
    yspan = y(N - n + 1:N);
    istart = j - N + n - 1;
    iend = j - N + n;
    yint(j) = yint(j-1)+localnc(tspan,yspan,istart,iend,n,p);
end
end

function y = localnc(tspan,yspan,istart,iend,n,p)
L = zeros(n,p);
X = zeros(n,1);
for i = 1:n
    %if i > 1
    %    X(i) = tspan(i) - tspan(i-1); %set x value proportional to time step
    %end
    X(i) = tspan(i) - tspan(1);
    for j = 1:p
        L(i,j) =  X(i)^(j-1);  %power of x, e.g. for 3/3: [1 X(1) X(1)^2;1 X(2) X(2)^2;1 X(3) X(3)^2];
    end
end

a = L\yspan; %solve SLAE or OLSM

xpowers  = zeros(1,p); % powers of x after integrating the polynomial, e.g. for P = 3 for [x x^2/2 x^3/3];
xmults = zeros(1,p);
for j = 1:p
    xpowers(j) = j;
    xmults(j) = 1/j;
end

%xstart = tspan(1);
%xend = tspan(iend); %set end point as you like

xstart = tspan(istart) - tspan(1);
xend = tspan(iend) - tspan(1); %set end point as you like

y = (xend.^xpowers.*xmults - xstart.^xpowers.*xmults) * a;
end
