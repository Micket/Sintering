% Creates a cube from (0,0,0) to (1,1,1)
close all
bbox = [-1,-1;1,1];
x = 1;
h = 0.13;
radius = 1.2;
a0 = sqrt(radius^2-x^2);
fixed_nodes = [ ...
                -x+a0,0;...
                 x-a0,0;...
                0, x-a0;...
                0,-x+a0;...
%                -x,0; x,0; 0,-x; 0,x;...
%                -x,-x; -x,x; x,-x; x,x;...
                ];

xx = [-2*x,-x,-x+a0,-x/2,0,x/2,x-a0,x,2*x]; % Same for yy
b = 1;
c = 2;
hh =[c,c,c,1,1,1,c,c,c;...
     c,c,c,1,1,1,c,c,c;...
     c,c,c,1,b,1,c,c,c;...
     1,1,1,1,1,1,1,1,1;...
     1,1,b,1,1,1,b,1,1;...
     1,1,1,1,1,1,1,1,1;...
     c,c,c,1,b,1,c,c,c;...
     c,c,c,1,1,1,c,c,c;...
     c,c,c,1,1,1,c,c,c];

%[p,t] = distmeshnd(@dsphere,@huniform,0.2,bbox,fixed_nodes, 0.5,0.5,0.5,0.5);
%[p,t] = distmeshnd(@d_4particle_rve,@(p) h_4particle_rve(p,1,radius),0.03,bbox,fixed_nodes, x,radius);
%[p,t] = distmesh2d(@d_4particle_rve,@huniform,0.02,bbox,fixed_nodes, x,radius);
[p,t] = distmeshnd(@d_4particle_rve, @(p)hmatrix(p,xx,xx,0,hh), h, bbox, fixed_nodes, x, radius);
%[p,t] = distmeshnd(@d_4particle_rve, @huniform, h, bbox, fixed_nodes, x, radius);

d = d_4particle_rve(p,x,radius);

nodes_left   = p(:,1) < -x+h/10;
nodes_right  = p(:,1) > x-h/10;
nodes_top    = p(:,2) > x-h/10;
nodes_bottom = p(:,2) < -x+h/10;
nodes_free   = d > -h/10 & ~nodes_left & ~nodes_right & ~nodes_top & ~nodes_bottom;

simpplot(p,t), hold on
plot(p(nodes_free,1),p(nodes_free,2),'or');
plot(p(nodes_left,1),p(nodes_left,2),'xr');
plot(p(nodes_right,1),p(nodes_right,2),'xb');
plot(p(nodes_top,1),p(nodes_top,2),'xg');
plot(p(nodes_bottom,1),p(nodes_bottom,2),'xy');

% Invert t, c = nodes connected to elements numbers
c = zeros(size(p,1),10);
q = ones(size(p,1),1);
for i = 1:size(t,1)
    for j = 1:3
        k = t(i,j);
        c(k,q(k)) = i;
        q(k) = q(k) + 1;
    end
end

tmp = c(nodes_left,:);   elem_left   = unique(tmp(tmp>0));
tmp = c(nodes_right,:);  elem_right  = unique(tmp(tmp>0));
tmp = c(nodes_top,:);    elem_top    = unique(tmp(tmp>0));
tmp = c(nodes_bottom,:); elem_bottom = unique(tmp(tmp>0));
tmp = c(nodes_free,:);   elem_free   = unique(tmp(tmp>0));

area = 0;
for i = 1:size(t,1)
    v1 = p(t(i,3),:) - p(t(i,1),:);
    v2 = p(t(i,2),:) - p(t(i,1),:);
    area = area + abs(det([v1;v2]))/2;
end
density = area/4

% Sanity check;
b = zeros(size(p,1),1);
b(t(:)) = 1;

nrelem = size(t,1)


%{
save 4particle_0.84_mesh.mat p t nodes_free nodes_left nodes_right nodes_top nodes_bottom elem_free elem_left elem_right elem_top elem_bottom
save 4particle_0.88_mesh.mat p t nodes_free nodes_left nodes_right nodes_top nodes_bottom elem_free elem_left elem_right elem_top elem_bottom
save 4particle_0.92_mesh.mat p t nodes_free nodes_left nodes_right nodes_top nodes_bottom elem_free elem_left elem_right elem_top elem_bottom
save 4particle_0.95_mesh.mat p t nodes_free nodes_left nodes_right nodes_top nodes_bottom elem_free elem_left elem_right elem_top elem_bottom
save 4particle_super_coarse_mesh.mat p t nodes_free nodes_left nodes_right nodes_top nodes_bottom elem_free elem_left elem_right elem_top elem_bottom
%}

