for radius = linspace(0.52,1/sqrt(2),100)
for n = 1:7

innerradius = 0.4;

inside = [-1, -1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2];
outside = [-1, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0];
bc = [0,0,0,0,1,1,1,1,1,1,1,1];
material = [1,2,3,4,5,5,5,5,5,5,5,5];
%eltype = [51, 51, 54, 54, 56, 56, 56, 56, 56, 56, 56, 56];
eltype = {'Tr21Stokes','Tr21Stokes','Line2SurfaceTension','Line2SurfaceTension','Line2BoundaryElement','Line2BoundaryElement','Line2BoundaryElement','Line2BoundaryElement','Line2BoundaryElement','Line2BoundaryElement','Line2BoundaryElement','Line2BoundaryElement'};
a = acos(0.5/radius);
rho = (radius^2*(pi/4-a) + sqrt(radius^2-0.5^2)/2)*4;

segments = 8 + 9*4*n + 8*n^2;

fprintf('radius %e, rho = %.2f\n',radius, rho);
basename = sprintf('rve_%d_%.2f',n,rho);

fid = fopen([basename,'.in'],'w');
fprintf(fid, ['%s\n',...
    'Simply 2D RVE with surface tension\n',...
    'StokesFlow nsteps 2000 lstype 3 smtype 7 nmodules 1 deltaT 0.0005 nonlinform 2 rtolv 1.0\n',...
    'vtk tstep_all domain_all primvars 2 4 5 cellvars 1 46 stype 0\n',...
    'domain 2dIncompFlow\n',...
    'OutputManager tstep_all dofman_all element_all\n',...
    'ndofman 0 nelem 0 ncrosssect 12 nmat 5 nbc 1 nic 0 nltf 1 topology particletopology\n'], basename);
fprintf(fid, 'EmptyCS %d\n',1:12);
fprintf(fid, ['NewtonianFluid 1 d 1 mu 1\n',...
    'NewtonianFluid 2 d 1 mu 10\n',...
    'SurfaceTension 3 g 1\n',...
    'SurfaceTension 4\n',...
    'SurfaceTension 5\n',...
    ...%'PrescribedGradient 1 loadTimeFunction 1 gradient 2 2 {0 0.0; 0.0 0} cCoord 2 %f %f defaultDofs 2 7 8\n',...
    'MixedGradientPressure 1 loadTimeFunction 1 devgradient 2 2 {0 0.0; 0.0 0} pressure 0.0 cCoord 2 %f %f defaultDofs 2 7 8\n',...
    'ConstantFunction 1 f(t) 1.0\n'], n/2, n/2);
eltypestring = '';
for i = 1:length(eltype)
    eltypestring = sprintf('%s %s',eltypestring,eltype{i});
end
fprintf(fid, ['ParticleTopology nsd 2 baseresolution %d bboxa 2 %f %f bboxb 2 %f %f neighbors %d tubewidth %f nsegments %d',...
  ' regioninside %d %s regionoutside %d %s material %d %s bc %d %s elementtype %d %s\n'],...
    n*250,-0.5,-0.5,n+0.5,n+0.5,5,1.1,...
    segments, length(inside), int2str(inside), length(outside), int2str(outside), ...
    length(material), int2str(material), length(bc), int2str(bc), length(eltype), eltypestring);

c = 1;
for i = 0:n-1
for j = 0:n-1
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 3, 0+i, 0+j, radius, a*180/pi, 90-a*180/pi ); c = c + 1;
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 3, 0+i, 1+j, radius, -90+a*180/pi, -a*180/pi ); c = c + 1;
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 3, 1+i, 1+j, radius, -180+a*180/pi, -90-a*180/pi); c = c + 1;
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 3, 1+i, 0+j, radius, 90+a*180/pi, 180-a*180/pi); c = c + 1;

    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 4, 0+i, 0+j, innerradius, 0, 90); c = c + 1;
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 4, 0+i, 1+j, innerradius, -90, 0); c = c + 1;
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 4, 1+i, 1+j, innerradius, -180, -90); c = c + 1;
    fprintf(fid, 'Circle %d id %d center 2 %f %f radius %f start %f end %f\n', c, 4, 1+i, 0+j, innerradius, 90, 180); c = c + 1;
end
end

for i = 0:n-1
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c, 5, 0, innerradius+i, 0, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c, 6, innerradius+i, n, 1-innerradius+i, n); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c, 7, n, 1-innerradius+i, n, innerradius+i); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c, 8, 1-innerradius+i, 0, innerradius+i, 0); c = c + 1;

    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 5, 0, innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 5, 0, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 6, innerradius+i, n); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 6, 1-innerradius+i, n); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 7, n, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 7, n, innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 8, 1-innerradius+i, 0); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 8, innerradius+i, 0); c = c + 1;

    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c, 9, 0, i, 0, innerradius+i); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c, 9, 0, 1-innerradius+i, 0, 1+i); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c,10, 0, n, innerradius+i, n); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c,10, 1-innerradius+i, n, 1+i, n); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c,11, n, 1+i, n, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c,11, n, innerradius+i, n, i); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c,12, 1+i, 0, 1-innerradius+i, 0); c = c + 1;
    fprintf(fid, 'Line %2d id %d start 2 %f %f end 2 %f %f\n', c,12, innerradius+i, 0, i, 0); c = c + 1;

    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 9, 0, innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 9, 0, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c,10, innerradius+i, n); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c,10, 1-innerradius+i, n); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c,11, n, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c,11, n, innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c,12, 1-innerradius+i, 0); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c,12, innerradius+i, 0); c = c + 1;

    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, 0, innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, 0, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, innerradius+i, n); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, 1-innerradius+i, n); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, n, 1-innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, n, innerradius+i); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, 1-innerradius+i, 0); c = c + 1;
    fprintf(fid, 'CornerPoint %2d id %d coords 2 %f %f\n', c, 4, innerradius+i, 0); c = c + 1;
end

fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c, 9, 0, 0); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c, 9, 0, n); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c,10, 0, n); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c,10, n, n); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c,11, n, n); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c,11, n, 0); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c,12, n, 0); c = c + 1;
fprintf(fid, 'CornerPoint %d id %d coords 2 %f %f\n', c,12, 0, 0); c = c + 1;

fclose(fid);

end
end
