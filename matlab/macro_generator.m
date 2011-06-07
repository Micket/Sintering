% Create a mesh of some sort.
n = 7;
[xx,yy] = meshgrid(linspace(-1,1,n),linspace(0,1,ceil(n/2)));
xx = xx(:);
yy = yy(:);

coords = [xx,yy];
elem = delaunay(coords(:,1),coords(:,2));


% Density function
r_min = 0.83;
r_max = 0.90;
rho = @(x) r_min + (r_max-r_min)*(1-x(2))*x(1)^2;


available = 0.83:0.01:1.00;

for i = 1:size(elem,1)
    gp = mean(coords(elem(i,1:3),:),1); % Gauss point (center)
    density = rho(gp);
    [x,j] = min(abs(available-density));
    filename = sprintf('rve_1_%.3f',available(j));

    fprintf('Density = %f -> %s\n',density,filename);
end

ndofman = size(coords,1);
nelem = size(elem,1);


% Write oofem input file.

fid = fopen('macro_problem_nonhom.in','w');

fprintf(fid, ['macro_problem_nonhom.out\n',...
'The macroscale problem, using sintering.\n',...
'IncrLinearStatic nmodules 2 deltat 0.002 nsteps 1000\n',...
'vtk tstep_all domain_all primvars 1 1 vars 3 1 4 75\n',...
'gp 1 tstep_step 1 domain_all ncoords 2 vars 3 1 4 75\n',...
'domain 2DPlaneStress\n',...
'OutputManager tstep_all dofman_all element_all\n']);

fprintf(fid, 'ndofman %d nelem %d ncrosssect 1 nmat %d nbc 1 nic 0 nltf 1\n',ndofman,nelem,length(available));

for i = 1:ndofman
    bc = [0,0];
    fprintf(fid, 'node %d coords 2 %f %f bc 2 %d %d\n',i,coords(i,1),coords(i,2),bc(1),bc(2));
end


for i = 1:nelem
    gp = mean(coords(elem(i,1:3),:),1); % Gauss point (center)
    density = rho(gp);
    [x,j] = min(abs(available-density));
    fprintf(fid, 'trplanestress2d %d nodes 3 %d %d %d crossSect 1 mat %d\n', i, elem(i,1), elem(i,2), elem(i,3), j);
end

fprintf(fid, 'SimpleCS 1 thick 1.0 width 1.0\n');

for i = 1:length(available)
    fprintf(fid, 'FE2SinteringMaterial 1 d 0.0 inputfile "/home/mikael/workspace/OOFEM/rve/rve_%.2f.in"\n', available(i));
end

fprintf(fid, 'BoundaryCondition 1 LoadTimeFunction 1 prescribedvalue 0.0\n');
fprintf(fid, 'ConstantFunction 1 f(t) 1.0\n');

fclose(fid);
