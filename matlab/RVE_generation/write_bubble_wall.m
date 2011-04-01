clear all

m = 2;
lpx = linspace(0,1,m);
[px,py] = meshgrid(lpx, lpx);
p = [px(:), py(:)];
if m == 2
    t = [1,3,4;1,4,2];
else
    t = delaunay(p(:,1), p(:,2));
end

n = size(t,1);
edges = sort([t(:,1),t(:,2); t(:,2),t(:,3); t(:,3),t(:,1)],2);
et = [1:n,1:n,1:n]';
[tmp,I] = sort(edges(:,1)*1e6 + edges(:,2));
edges = edges(I,:);
et = et(I);

% Remove duplicates
et = [et,zeros(size(et))];
remove = zeros(size(edges,1),1);
for i = 2:size(edges,1)
    if edges(i,:) == edges(i-1,:)
        remove(i) = true;
        et(i-1,2) = et(i,1);
    end
end
%remove(I) = remove;
%edges(I,:) = edges;
%et(I,:) = et;
edges = edges(~remove,:);
et = et(~remove,:);

% Invert triangle->element list (et)
te = zeros(n,3);
k = zeros(n,1);
for i = 1:size(et,1)
    for j = 1:2
        q = et(i,j);
        if q > 0
            k(q) = k(q) + 1;
            te(q,k(q)) = i;
        end
    end
end

% Fix ordering... sigh, there must be some nicer way than all this shit.
for i = 1:size(t,1)
    ee = edges(te(i,:),:); % already sorted
    tt1 = sort(t(i,[1,2]));
    tt2 = sort(t(i,[2,3]));
    tt3 = sort(t(i,[3,1]));
    order = [find(ee(:,1) == tt1(1) & ee(:,2) == tt1(2)),...
             find(ee(:,1) == tt2(1) & ee(:,2) == tt2(2)),...
             find(ee(:,1) == tt3(1) & ee(:,2) == tt3(2))];
    te(i,:) = te(i,order);
end

% Add nodes
new_p = (p(edges(:,1),:) + p(edges(:,2),:))/2;

% Construct new elements
p2 = [p;new_p];
t2 = [t,te+size(p,1)];
edges2 = [edges,(1:size(edges,1))'+size(p,1)];

new_p_left   = new_p(:,1) < min(new_p(:,1))+1e-5;
new_p_right  = new_p(:,1) > max(new_p(:,1))-1e-5;
new_p_bottom = new_p(:,2) < min(new_p(:,2))+1e-5;
new_p_top    = new_p(:,2) > max(new_p(:,2))-1e-5;

p2_left   = p2(:,1) < min(p2(:,1))+1e-5;
p2_right  = p2(:,1) > max(p2(:,1))-1e-5;
p2_bottom = p2(:,2) < min(p2(:,2))+1e-5;
p2_top    = p2(:,2) > max(p2(:,2))-1e-5;

edges_left = find(new_p_left);
edges_bottom = find(new_p_bottom);
edges_top = find(new_p_top);
edges_right = find(new_p_right);

edges_free = find(new_p_right | new_p_top);

if 0
clf
hold on

for i = 1:size(t2,1)
    plot(p2(t2(i,[1,4,2,5,3,6,1]),1),p2(t2(i,[1,4,2,5,3,6,1]),2),'k');
end
for i = 1:size(edges_left,1)
    plot(p(edges(edges_left(i),:),1),p(edges(edges_left(i),:),2),'b','linewidth',2);
end
for i = 1:size(edges_bottom,1)
    plot(p(edges(edges_bottom(i),:),1),p(edges(edges_bottom(i),:),2),'r','linewidth',2);
end
for i = 1:size(edges_free,1)
    plot(p(edges(edges_free(i),:),1),p(edges(edges_free(i),:),2),'m');
end
end

explanation = 'Simply 2D drop resting in a corner, driven by surface tension';
fname = 'bubble_wall';
gamma_s = [0.2,0.05,0.0];
nsteps = 100;
dt = 1;

theta = acos((gamma_s(3) - gamma_s(2))/gamma_s(1))*180/pi


ndofman = size(p2,1);
nelem = 0;
nodes = [];

for i = 1:size(p2,1)
    if i <= size(p,1)
        dofidmask = [7,8,11];
        BC = [0,0,0];
    else
        dofidmask = [7,8];
        BC = [0,0];
    end
    if p2_left(i) 
        BC(1) = 1;
    end
    if p2_bottom(i)
        BC(2) = 1;
    end

    nodes = [nodes,...
            sprintf('node %4d coords 2 %17.10e %17.10e ndofs %d DofIDMask %d %s bc %d %s\n',...
            i,p2(i,:),size(dofidmask,2),size(dofidmask,2),num2str(dofidmask),size(BC,2),num2str(BC))];
end

nelem = 0;
elements = sprintf('Tr21Stokes %4d nodes 6 %4d %4d %4d %4d %4d %4d crossSect 1 mat 1\n',[(1:size(t2,1))',t2]');
nelem = nelem + size(t2,1);

elements = [elements,...
    sprintf('Line2SurfaceTension %3d nodes 3 %4d %4d %4d crossSect 2 mat 2\n',[(1:size(edges_top,1))'+nelem,edges2(edges_top,:)]')];
nelem = nelem + size(edges_top,1);

elements = [elements,...
    sprintf('Line2SurfaceTension %3d nodes 3 %4d %4d %4d crossSect 2 mat 2\n',[(1:size(edges_right,1))'+nelem,edges2(edges_right,:)]')];
nelem = nelem + size(edges_right,1);

elements = [elements,...
    sprintf('Line2SurfaceTension %3d nodes 3 %4d %4d %4d crossSect 3 mat 3\n',[(1:size(edges_bottom,1))'+nelem,edges2(edges_bottom,:)]')];
nelem = nelem + size(edges_bottom,1);

elements = [elements,...
    sprintf('Line2SurfaceTension %3d nodes 3 %4d %4d %4d crossSect 5 mat 5\n',[(1:size(edges_left,1))'+nelem,edges2(edges_left,:)]')];
nelem = nelem + size(edges_left,1);

material = sprintf(['NewtonianFluid 1 d 1.0 mu 1.0\n'...
                    'SurfaceTension 2 g %e\n',...
                    'SurfaceTension 3 g %e\n'...
                    'SurfaceTension 4 g %e\n'...
                    'SurfaceTension 5 \n'],gamma_s);

boundarycondition = sprintf('BoundaryCondition 1 LoadTimeFunction 1 PrescribedValue 0.0\n');

loadtimefunc = sprintf('ConstantFunction 1 f(t) 1.0\n');

crosssections = sprintf('EmptyCS %d\n',1:5);

header = sprintf([...
    '%s.out\n%s\n'...
    'StokesFlow nsteps %d lstype 3 smtype 7 nmodules 1 deltaT %f nonlinform 2 rtolv %e\n',...
    'vtk tstep_all domain_all primvars 2 4 5 cellvars 1 46\n',...
    'domain 2dIncompFlow\n',...
    'OutputManager tstep_all dofman_all element_all\n',...
    'ndofman %d nelem %d ncrosssect 5 nmat 5 nbc 1 nic 0 nltf 1\n'],...
    fname, explanation, nsteps, dt, 1, ndofman, nelem);

fid = fopen([fname,'.in'],'w+');
fwrite(fid,header);
fwrite(fid,nodes);
fwrite(fid,elements);
fwrite(fid,crosssections);
fwrite(fid,material);
fwrite(fid,boundarycondition);
fwrite(fid,loadtimefunc);
system(['cp ',fname,'.in ~/workspace/OOFEM/']);

