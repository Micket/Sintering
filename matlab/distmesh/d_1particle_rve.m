function d = d_1particle_rve(p,x,r)

a0 = sqrt(r^2-x^2);

ds = drectangle0(p,0,x,0,x);
dc = dcircle(p,x,x,r);

d = max(dc,ds); % Intersection
