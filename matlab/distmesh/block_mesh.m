clear all
close all

% Creates a cube from (0,0,0) to (1,1,1)
bbox = [-1,-1;1,1];
%bbox = [-1,-1,-1;2,2,2];

fixed_nodes = [];

%[p,t] = distmeshnd(@dsphere,@huniform,0.2,bbox,fixed_nodes, 0.5,0.5,0.5,0.5);
%[p,t] = distmeshnd(@dblock0,@huniform,0.2,bbox,fixed_nodes, 0,1,0,1,0,1);

[p,t] = distmesh2d(@dcircle_special,@huniform,0.1,bbox,fixed_nodes, 0,0,1);

simpplot(p,t,'p(:,2)>0');

