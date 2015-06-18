% Case D_3:
m = 3;
hm = 0.57735; % Obtained by sqrt(tan(maximum_angle/2)^2 - 1)/2) (I didn't understand it's definition)

C = @(a,b,c) acos((cos(c) - cos(a) .* cos(b))./(sin(a).*sin(b)));
S2 = @(a, b, c) 2*(pi - C(a, b, c) - cos(a) .* C(c, a, b) - cos(b) .* C(b, c, a));
S1 = @(a) 2 * pi * (1 - cos(a));

a1 = @(r) acos(hm./r);
a2 = @(r) acos(1./r);


h = @(x) x > 0;

x = @(r) 4*pi ...
	- h(r-hm) .* (2 *S1(a1(r)))  ...
	- h(r-1) .* (2*m*S1(a2(r))) ...
	+ h(r-sqrt(1+hm^2)) .* (4*m*S2(a1(r), a2(r), pi/2) + 2*m*S2(a2(r), a2(r), pi/m));

K = 1/9.5488e+00; % Don't know what this should be either.
p = @(w) K / (2*pi^2) .* sin(w/2).^2 .* x(tan(w/2));


ws = linspace(0.01, 104 * pi / 180, 100);
plot(ws*180/pi, p(ws));

% trapz(ws*180/pi, p(ws)) Ensure that the value is normed. (I don't know how to compute K analytically)

data = [ws'*180/pi, real(p(ws))'];
save('D3.txt', '-ascii', 'data')
