n = 5;

x = 0:n;

[xx,yy] = meshgrid(x,x);
xx = xx(:);
yy = yy(:);

tri = delaunay(xx,yy);

cx = mean(xx(tri),2); cy = mean(yy(tri),2);
sortvalue = cy*1000 + cx;
[crap,index] = sort(sortvalue);

tri = tri(index,:);

fprintf('ndofman %d nelem %d ncrosssect 1 nmat 1 nbc 1 nic 0 nltf 1\n',length(xx),size(tri,1));
fprintf('node %d coords 2 %d %d\n',[(1:length(xx))',xx,yy]');
fprintf('TrPlaneStress %2d nodes 3 %2d %2d %2d crossSect 1 mat 1\n',[(1:size(tri,1))',tri]');

cx = mean(xx(tri),2); cy = mean(yy(tri),2);


%t = linspace(0,3,301); fprintf('%.2f ',t);
