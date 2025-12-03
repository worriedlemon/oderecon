function h_ji = h_ji_byint(Ej,t,xi)

Es = integrate_trapz(Ej,0,t);
Es0 = trapz(t,Es)/(t(end) - t(1)); % calculate bias
Es = Es - Es0;

h_ji =  trapz(t, Es' .* xi);