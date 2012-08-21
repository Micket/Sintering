clear
% Create a mesh of some sort.
n = 3;
m = n-1; % Number of partitions (parallel computations)
cells = 1;
available = 0.83:0.01:1.00;

deltaT = 0.003;
nsteps = 5;

[xx,yy] = meshgrid(linspace(-1,1,n),linspace(0,1,ceil(n/2)));
coords = [xx(:),yy(:)];
elem = delaunay(coords(:,1),coords(:,2));
[elem, coords] = lin2quad(elem,coords);
xx = coords(:,1); yy = coords(:,2);
cx = mean(xx(elem),2);

nodepartitions = cell(size(coords,1),1);
numnodes = zeros(m,1);

for i = 1:m
    limit = -1.01+2.02*i/m;
    oldlimit = -1.01+2.02*(i-1)/m;
    pelem = find(cx<=limit & cx > oldlimit);
    partitions{i} = pelem;
    for j = 1:length(pelem)
        for k = 1:length(elem(pelem,:))
            q = nodepartitions{elem(pelem(j),k)};
            if length(intersect(q,i)) == 0
                nodepartitions{elem(pelem(j),k)} = [q,i];
                numnodes(i) = numnodes(i) + 1;
            end
        end
    end
end

% Density function
r_min = 0.83;
r_max = 0.97;
rho = @(x) r_min + (r_max-r_min)*(1-x(2))*(0+x(1)^2);

ndofman = size(coords,1);
nelem = size(elem,1);


% Write oofem input file.
for mi = 1:m
pelem = partitions{mi};

fid = fopen(['macro_problem_nonhom.in.',num2str(mi-1)],'w');

fprintf(fid, ['macro_problem_nonhom\n',...
'The macroscale problem, using sintering.\n',...
'StokesFlow nmodules 1 deltat %e nsteps %d nonlinform 2 smtype 7 lstype 3 rtolv 1e-6 rtold 1e6 manrmsteps 1\n',...
'vtkxml tstep_all domain_all primvars 2 4 5 cellvars 2 43 76\n',...
'domain 2DPlaneStress\n',...
'OutputManager tstep_all dofman_all element_all\n'], deltaT, nsteps);

fprintf(fid, 'ndofman %d nelem %d ncrosssect 1 nmat %d nbc 1 nic 0 nltf 1\n',numnodes(mi),length(pelem),length(available));

for i = 1:ndofman
    if i <= n*ceil(n/2)
        dofidmask = [7,8,11];
        bc = [0,0,0];
    else 
        dofidmask = [7,8];
        bc = [0,0];
    end

    if abs(coords(i,1)) < eps
        %if coords(i,2) < 1/(ceil(n/2))*2 + eps
        bc(1) = 1;
        %end
    end
    if coords(i,2) < eps
        bc(2) = 1;
    end

    np = nodepartitions{i};
    if ~ismember(mi,np)
        continue
    end
    if length(np) > 1
        fprintf(fid, 'node %2d coords 2 %9.6f %9.6f', i, coords(i,1), coords(i,2));
        fprintf(fid, ' bc %d %s',length(bc), num2str(bc));
        fprintf(fid, ' ndofs %d dofidmask %d %s',length(dofidmask), length(dofidmask), num2str(dofidmask));
        fprintf(fid, ' shared partitions %d %s\n', length(np), num2str(np-1));
    else
        fprintf(fid, 'node %2d coords 2 %9.6f %9.6f', i, coords(i,1), coords(i,2));
        fprintf(fid, ' bc %d %s',length(bc), num2str(bc));
        fprintf(fid, ' ndofs %d dofidmask %d %s\n',length(dofidmask), length(dofidmask), num2str(dofidmask));
    end
end


for i = 1:length(pelem)
    k = pelem(i);
    gp = mean(coords(elem(k,1:3),:),1); % Gauss point (center)
    density = rho(gp);
    [x,j] = min(abs(available-density));
    fprintf(fid, 'tr21stokes %d nodes %d ', k, length(elem(k,:)));
    fprintf(fid, '%d ', elem(k,:));
    fprintf(fid, 'crossSect 1 mat %d\n', j);
end
fprintf(fid, 'SimpleCS 1 thick 1.0 width 1.0\n');

for i = 1:length(available)
    %fprintf(fid, 'FE2SinteringMaterial %d d 0.0 inputfile "/beda/users/home/ohmanm/rve/rve_%d_%.2f.in"\n', i, cells, available(i));
    fprintf(fid, 'FE2FluidMaterial %d d 0.0 inputfile "/home/mikael/rve/rve_%d_%.2f.in"\n', i, cells, available(i));
    %fprintf(fid, 'NewtonianFluid %d d 0.0 mu 1.0\n', i);
end

fprintf(fid, 'BoundaryCondition 1 LoadTimeFunction 1 prescribedvalue 0.0\n');
fprintf(fid, 'ConstantFunction 1 f(t) 1.0\n');

fclose(fid);

end
