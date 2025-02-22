function s = intdiff4(x,th)
% -- s = intdiff4(xs,ths)
% Returnes integral over evenly sampled in time
% s = int(th*x'*dt),
% using 4 order numerical integration method and asymmetric differentiation to handle possible noise in data
% input:
% xs  is an N x 1 vector
% ths is an N x 1 vector
% output:
% s is a scalar


n = length(x);
s = 0;
for i = 2:n
    if i <= 3

        dx3 =    -x(i-1)/3 +  3*x(i)/2  - 3*x(i+1) + 11*x(i+2)/6;
        dx2 =     x(i-1)/6    - x(i)      + x(i+1)/2  + x(i+2)/3;
        dx1 =    -x(i-1)/3    - x(i)/2    + x(i+1)    - x(i+2)/6;
        dx0 = -11*x(i-1)/6  + 3*x(i)    - 3*x(i+1)/2  + x(i+2)/3;

        f3 = dx3 * th(i + 2);
        f2 = dx2 * th(i + 1);
        f1 = dx1 * th(i);
        f0 = dx0 * th(i - 1);
        
        s = s + (3*f0)/8 + (19*f1)/24 - (5*f2)/24 + f3/24;

    elseif i >= n - 1

        dx3 =    -x(i-3)/3 +  3*x(i-2)/2  - 3*x(i-1) + 11*x(i)/6;
        dx2 =     x(i-3)/6    - x(i-2)      + x(i-1)/2  + x(i)/3;
        dx1 =    -x(i-3)/3    - x(i-2)/2    + x(i-1)    - x(i)/6;
        dx0 = -11*x(i-3)/6  + 3*x(i-2)    - 3*x(i-1)/2  + x(i)/3;

        f3 = dx3 * th(i);
        f2 = dx2 * th(i - 1);
        f1 = dx1 * th(i - 2);
        f0 = dx0 * th(i - 3);

        s = s + f0/24 - (5*f1)/24 + (19*f2)/24 + (3*f3)/8;
    else

        dx3 =    -x(i-2)/3 +  3*x(i-1)/2  - 3*x(i) + 11*x(i+1)/6;
        dx2 =     x(i-2)/6    - x(i-1)      + x(i)/2  + x(i+1)/3;
        dx1 =    -x(i-2)/3    - x(i-1)/2    + x(i)    - x(i+1)/6;
        dx0 = -11*x(i-2)/6  + 3*x(i-1)    - 3*x(i)/2  + x(i+1)/3;

        f3 = dx3 * th(i+1);
        f2 = dx2 * th(i);
        f1 = dx1 * th(i - 1);
        f0 = dx0 * th(i - 2);

        s = s + (-1/24*f3 + 13/24*f2 + 13/24*f1 - 1/24*f0);
    end
end

end