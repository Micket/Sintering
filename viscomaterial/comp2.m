clear all 
close all

% load turbine p e t from pde-tool, can be downloaded from Study portal
%load turbine_init
load turbine_fine
% Create Coord, Dof and Edof from p, t 
nnod = max(size(p));
nrelem = max(size(t));
nrdof = nnod*2;
Coord=p';
Dof = reshape(1:nrdof,2,nnod)';
Edof = [(1:nrelem)',Dof(t(1,:),:),Dof(t(2,:),:),Dof(t(3,:),:)];
[Ex,Ey] = coordxtr(Edof,Coord,Dof,3);

% set x displ of these nodes and y displ for the first node to = 0
nod_left = find(Coord(:,1) < 1e-6);
bc = [Dof(nod_left(:),1), zeros(size(nod_left)); Dof(nod_left(1),2), 0];
[dummy,outmost] = max(Coord(:,1)); % Find a node on the right edge

t = 6;    % plate thickness (element thickness)
gam = [eye(3),zeros(3);zeros(3),eye(3)/2;zeros(3),eye(3)/2];

% loads and other constants
tmax      = 1000;
ntsteps   = 100; % total number of timesteps
time_init = 100; % initialization of load via ramp until time_init
dtmax     = tmax/ntsteps;
dtmin     = 0.1*tmax/ntsteps;

errtol    = 1e-5; % tolerance for equilibrium iterations
ptype     = 2;    % plane strain
magn      = 75;   % magnification factor for deformed mesh
maxiter   = 4;   % max number of iterations
r_0       = 800;  % volume force 
force_0   = 8.2/r_0;

modpar.E  = 150e3; 
modpar.H  = 150e3;
modpar.nu = 0.3;
modpar.sy = 10e0;
%modpar.sy = 0e0;
modpar.sy_dyn = 1000e0;
modpar.n = 4;
modpar.S = 1.5e-2;
modpar.m = 2;
modpar.k_inf = 3e2;
modpar.kappa_inf_h = 3e2;
modpar.kap_inf = 3e2;
modpar.kappa_inf = 3e2;
modpar.t_star = 1e4;
modpar
De = hooke(4,modpar.E,modpar.nu);

% Initiation of variables
u         = zeros(nrdof,1);
u_hist    = zeros(ntsteps,1);
t_hist    = zeros(ntsteps,1);
dam_hist  = zeros(ntsteps,1); 
sig       = zeros(9,nrelem);
eps       = zeros(9,nrelem);
oldsig    = zeros(9,nrelem); 
oldeps    = zeros(9,nrelem);
lconv     = ones(nrelem,1);
%state_old = struct('eps_e',zeros(9,1),'k',0,'b',1,'mu_h',num2cell(zeros(nrelem,1)));
state_old = struct('eps_e',zeros(9,1),'k',0,'b',1,'mu',num2cell(zeros(nrelem,1)));
state_new = state_old;

%% Global timestepping
stop_iteration = false;
time = 0;
dt = dtmax;
du = zeros(nrdof,1);
for tstep = 1:ntsteps
    iter = 0;
    olderr = inf;
    for iter = 1:maxiter
        Ded = extract(Edof,du);
        % This part is does the same thing as assem, only *much* faster.
        apos = 0; II = zeros(nrelem*6^2,1); JJ = zeros(nrelem*6^2,1); XX = zeros(nrelem*6^2,1);
        fint = zeros(nrdof,1); 
        fext = zeros(nrdof,1);
        for i=1:nrelem
            [dummy,Deps_cf]=plants(Ex(i,:),Ey(i,:),[ptype t],De,Ded(i,:)); 
            Deps = gam*Deps_cf';

            [eps(:,i),sig(:,i),state_new(i),E_tang,lconv(i)] = ...
              visco(oldeps(:,i),oldsig(:,i),state_old(i),Deps,dt,modpar);

            sig_cf      = sig(1:6,i)';
            E_tang_cf   = E_tang(1:6,1:6); E_tang_cf(4:6,4:6) = E_tang_cf(4:6,4:6)/2;
            radius      = sum(Ex(i,:))/3+r_0;
            bforce      = [force_0*radius; 0]*min(1,(time+dt)/time_init);
            [Ke,fext_e] = plante(Ex(i,:),Ey(i,:),[ptype t],E_tang_cf,bforce);
            fint_e      = plantf(Ex(i,:),Ey(i,:),[ptype t],sig_cf)'-fext_e;

            % The loops do the same as assem
            for ii = 1:6, for jj = 1:6
            apos = apos+1; 
            II(apos) = Edof(i,1+ii); JJ(apos) = Edof(i,1+jj); XX(apos) = Ke(ii,jj);
            end, end
            fint(Edof(i,2:end)) = fint(Edof(i,2:end)) + fint_e;
            fext(Edof(i,2:end)) = fext(Edof(i,2:end)) + fext_e;
        end
        Kglob = sparse(II,JJ,XX); % See above. This is how sparse should be used.
        Q = zeros(nrdof,1); Q(bc(:,1)) = fint(bc(:,1));
        res = Q-fint;
        err = norm(res)/norm(fext);

        if err > olderr
            fprintf('Err not decreasing!\n')
            dt = dt/2;
            du = (du-ddu)/2;
            if dt < dtmin
                stop_iteration = true;
            end
            olderr = inf;
            continue;
        %elseif err < olderr
        %    dt = dt*1.5;
        %    du = (du-ddu)*1.5;
        end
        olderr = err;

        ddu = solveq(Kglob,res,bc);
        du = du+ddu;
        fprintf('tstep %3d iter %2d, dt = %.2f, err = %.4e b = %.2f\n',...
            tstep,iter,dt,err,min([state_new.b]))
        if err < errtol
            break
        end
    end
    if ~all(lconv) || stop_iteration
        break;
    end
    u=u+du;
    time = time + dt;
    oldsig = sig;
    oldeps = eps;
    state_old = state_new;
    u_hist(tstep)   = u(Dof(outmost,1));
    t_hist(tstep)   = time;
    dam_hist(tstep) = 1-min([state_new.b]);
end 



% Plotting
ed = extract(Edof,u);
figure, hold on, axis equal
eldraw2(Ex,Ey,[1 3 0]);
eldisp2(Ex,Ey,ed,[1 1 0],magn); 
xlabel('x'), ylabel('y')

figure, hold on, axis equal, colorbar
patch(Ex',Ey',1-[state_new.b])
title('Damage','interpreter','latex','fontsize',14)

% Skip this in your report
figure, hold on, axis equal, colorbar
patch(Ex',Ey',1-[state_new.k])
title('Internal variable $k$','interpreter','latex','fontsize',14)

% Skip this figure in your report
sig_dev = sig; sig_dev(1:3,:) = sig(1:3,:) - [1;1;1]*mean(sig(1:3,:));
sig_vm = sqrt(sum(sig_dev.*sig_dev)*3/2);
figure, hold on, axis equal, colorbar
patch(Ex',Ey',sig_vm)
title('Effective stress','interpreter','latex','fontsize',14)

figure, grid on
plot(t_hist(1:tstep-1),u_hist(1:tstep-1))
ylabel('Displacement [mm]','interpreter','latex','fontsize',18)
xlabel('$t$ [s]','interpreter','latex','fontsize',18)

figure, grid on
plot(t_hist(1:tstep-1),dam_hist(1:tstep-1))
ylabel('Damage','interpreter','latex','fontsize',18)
xlabel('$t$ [s]','interpreter','latex','fontsize',18)

