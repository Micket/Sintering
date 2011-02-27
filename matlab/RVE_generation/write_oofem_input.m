clear all
close all

% 0.84, 0.88, 0.92, 0.95
fname = '4particle_0.84_mesh';
%fname = '4particle_0.88_mesh';
%fname = '4particle_0.92_mesh';
%fname = '4particle_0.95_mesh';
%fname = '4particle_coarse_mesh';
%fname = '4particle_super_coarse_mesh';
q1_to_q2

explanation = 'Simply 2D RVE with surface tension';
gamma_s = 0.2;
ndofman = size(p2,1);
nelem = 0;
r = 0.82; % radius of the embedded, stiff particle

nodes = [];
boundary = zeros(size(p2,1),1);
boundary(edges2(edges_free,:)) = 1;
boundary(bc_nodes2) = 1;

for i = 1:size(p2,1)
    if i <= size(p,1)
        dofidmask = [7,8,11];
        BC = [0,0,0];
    else
        dofidmask = [7,8];
        BC = [0,0];
    end
    %if nodes_left2(i);   BC(1) = 1; end
    %if nodes_right2(i);  BC(2) = 2; end
    %if nodes_top2(i);    BC(1) = 3; end
    %if nodes_bottom2(i); BC(2) = 4; end
    if bc_nodes2(i); BC(1:2) = 1; end

    bflag = '';
    if boundary(i)
        bflag = 'boundary';
    end

    nodes = [nodes,...
            sprintf('node %4d coords 2 %17.10e %17.10e ndofs %d DofIDMask %d %s bc %d %s %s\n',...
            i,p2(i,:),size(dofidmask,2),size(dofidmask,2),num2str(dofidmask),size(BC,2),num2str(BC),bflag)];
end

q2q1trstokes = '';
surfacetension = '';

px = p(:,1); py = p(:,2);
cx = mean(px(t),2); cy = mean(py(t),2);
x1 = px(t(:,1)); x2 = px(t(:,2)); x3 = px(t(:,3));
y1 = py(t(:,1)); y2 = py(t(:,2)); y3 = py(t(:,3));

is_particle = sqrt((x1-1).^2+(y1-1).^2) < r | sqrt((x1+1).^2+(y1-1).^2) < r | sqrt((x1-1).^2+(y1+1).^2) < r | sqrt((x1+1).^2+(y1+1).^2) < r | ...
              sqrt((x2-1).^2+(y2-1).^2) < r | sqrt((x2+1).^2+(y2-1).^2) < r | sqrt((x2-1).^2+(y2+1).^2) < r | sqrt((x2+1).^2+(y2+1).^2) < r | ...
              sqrt((x3-1).^2+(y3-1).^2) < r | sqrt((x3+1).^2+(y3-1).^2) < r | sqrt((x3-1).^2+(y3+1).^2) < r | sqrt((x3+1).^2+(y3+1).^2) < r;
mat = 1+is_particle;
patch([x1,x2,x3]',[y1,y2,y3]',mat')
for i = 1:size(t2,1)
    q2q1trstokes = [q2q1trstokes,...
        sprintf('Tr21Stokes %4d nodes 6 %3d %3d %3d %3d %3d %3d crossSect 1 mat %d\n',i+nelem,t2(i,:),mat(i))];
end
nelem = nelem + size(t2,1);

if size(edges_free,1) > 0
    surfacetension = sprintf('Line2SurfaceTension %3d nodes 3 %3d %3d %3d crossSect 2 mat 3\n',...
        [(1:size(edges_free,1))'+nelem,edges2(edges_free,:)]');
end
nelem = nelem + size(edges_free,1);

% Edges for topology description;
if size(edges_left,1) > 0
    surface_left = sprintf('Line2SurfaceTension %3d nodes 3 %3d %3d %3d crossSect 3 mat 4\n',...
        [(1:size(edges_left,1))'+nelem,edges2(edges_left,:)]');
end
nelem = nelem + size(edges_left,1);

% Edges for topology description;
if size(edges_top,1) > 0
    surface_top = sprintf('Line2SurfaceTension %3d nodes 3 %3d %3d %3d crossSect 4 mat 4\n',...
        [(1:size(edges_top,1))'+nelem,edges2(edges_top,:)]');
end
nelem = nelem + size(edges_top,1);

% Edges for topology description;
if size(edges_right,1) > 0
    surface_right = sprintf('Line2SurfaceTension %3d nodes 3 %3d %3d %3d crossSect 5 mat 4\n',...
        [(1:size(edges_right,1))'+nelem,edges2(edges_right,:)]');
end
nelem = nelem + size(edges_right,1);

% Edges for topology description;
if size(edges_bottom,1) > 0
    surface_bottom = sprintf('Line2SurfaceTension %3d nodes 3 %3d %3d %3d crossSect 6 mat 4\n',...
        [(1:size(edges_bottom,1))'+nelem,edges2(edges_bottom,:)]');
end
nelem = nelem + size(edges_bottom,1);

material = sprintf(['NewtonianFluid 1 d 1.0 mu 1.0\n'...
                    'NewtonianFluid 2 d 1.0 mu 10.0\n',...
                    'SurfaceTension 3 g %e\n',...
                    'SurfaceTension 4\n'],gamma_s);

grad = [0,0.1;0.1,0];
dt = 0.1;
rtolv = 1;
nsteps = 100;

boundarycondition = sprintf('PrescribedGradient 1 loadTimeFunction 1 cCoord 2 0.0 0.0 gradient 2 2 {%e %e; %e %e}\n',grad');

loadtimefunc = sprintf('ConstantFunction 1 f(t) 1.0\n');


header = sprintf([...
    '%s.out\n%s\n'...
    'StokesFlow nsteps %d lstype 3 smtype 7 nmodules 1 deltaT %f nonlinform 2 rtolv %e\n',...
    'vtk tstep_all domain_all primvars 2 4 5 cellvars 1 46\n',...
    'domain 2dIncompFlow\n',...
    'OutputManager tstep_all dofman_all element_all\n',...
    'ndofman %d nelem %d ncrosssect 6 nmat 4 nbc %d nic 0 nltf 1\n'],...
    fname, explanation, nsteps, dt, rtolv, ndofman, nelem, 1);
 
fid = fopen([fname,'.in'],'w+');

fwrite(fid,header);
fwrite(fid,nodes);
fwrite(fid,q2q1trstokes);
fwrite(fid,surfacetension);
fwrite(fid,surface_left);
fwrite(fid,surface_top);
fwrite(fid,surface_right);
fwrite(fid,surface_bottom);
fwrite(fid,sprintf('EmptyCS 1\nEmptyCS 2\nEmptyCS 3\nEmptyCS 4\nEmptyCS 5\nEmptyCS 6\n'));
fwrite(fid,material);
fwrite(fid,boundarycondition);
fwrite(fid,loadtimefunc);

system(['cp ',fname,'.in ~/workspace/OOFEM/']);

