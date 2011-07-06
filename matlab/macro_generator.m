% Create a mesh of some sort.
n = 7;
[xx,yy] = meshgrid(linspace(-1,1,n),linspace(0,1,ceil(n/2)));
xx = xx(:);
yy = yy(:);
coords = [xx,yy];
elem = delaunay(coords(:,1),coords(:,2));


% Density function
r_min = 0.83;
r_max = 0.96;
rho = @(x) r_min + (r_max-r_min)*(1-x(2))*(0+x(1)^2);

cells = 2;
available = 0.83:0.01:1.00;
ndofman = size(coords,1);
nelem = size(elem,1);


% Write oofem input file.
fid = fopen('macro_problem_nonhom.in','w');

fprintf(fid, ['macro_problem_nonhom\n',...
'The macroscale problem, using sintering.\n',...
'IncrLinearStatic nmodules 1 deltat 0.002 nsteps 1000\n',...
'vtk tstep_all domain_all primvars 1 1 cellvars 1 43 stype 0\n',...
'domain 2DPlaneStress\n',...
'OutputManager tstep_all dofman_all element_all\n']);

fprintf(fid, 'ndofman %d nelem %d ncrosssect 1 nmat %d nbc 1 nic 0 nltf 1\n',ndofman,nelem,length(available));

for i = 1:ndofman
    bc = [0,0];
    if abs(coords(i,1)) < eps
        if coords(i,2) < 1/(ceil(n/2))*2 + eps
           bc(1) = 1;
        end
        if coords(i,2) < eps
            bc(2) = 1;
        end
    end
    fprintf(fid, 'node %d coords 2 %f %f bc 2 %d %d\n',i,coords(i,1),coords(i,2),bc(1),bc(2));
end


for i = 1:nelem
    gp = mean(coords(elem(i,1:3),:),1); % Gauss point (center)
    density = rho(gp);
    [x,j] = min(abs(available-density));
    disp(available(j))
    fprintf(fid, 'trplanestress2d %d nodes 3 %d %d %d crossSect 1 mat %d\n', i, elem(i,1), elem(i,2), elem(i,3), j);
end

fprintf(fid, 'SimpleCS 1 thick 1.0 width 1.0\n');

for i = 1:length(available)
    fprintf(fid, 'FE2SinteringMaterial %d d 0.0 inputfile "/home/mikael/workspace/OOFEM/rve/rve_%d_%.2f.in"\n', i, cells, available(i));
end

fprintf(fid, 'BoundaryCondition 1 LoadTimeFunction 1 prescribedvalue 0.0\n');
fprintf(fid, 'ConstantFunction 1 f(t) 1.0\n');

fclose(fid);
