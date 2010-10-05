clear all
close all

% Parameters
maxiter = 10;
n = 1000; % should just be large enough. 
tol = 1e-6; % Internal forces relative the reaction forces.
d  = 2; % Number of dimensions to run
E  = 200e9; 
nu = 0.3;
tau_y = 400e6;
tau_c = tau_y; % TODO
H = E/20;
t_star = 1;
eta0 = 1;
Dt = 1;
n = 2;
c = tau_y*0.5;

mu = E/(2*(1+nu));
lambda = E*nu/((1-2*nu)*(1+nu));
modpar = struct('mu',mu,'lambda',lambda,'tau_y',tau_y,'tau_c',tau_c,'H',H,'t_star',t_star,'eta0',eta0,'n',n,'c',c);

t = 1e-3;
%umax = 1e-4;
umax = 5e-3;
dux0 = 1e-6; % Just something small should suffice
dux = dux0;
max_dux = 1e-3;
min_dux = 1e-7;


% Preperation
load cass2_mesh_data.mat
%load small_mesh_data.mat
Ex = Ex*1e-3; Ey = Ey*1e-3;
Coord = Coord*1e-3;
dof_prescr = dof_prescr'; dof_fixed = dof_fixed'; dof_free = dof_free';
miny = min(Coord(:,2));
dof_vert = Dof(Coord(:,2)<=miny+10*eps(miny),1);
nrdof = numel(Dof);
nrelem = size(Edof,1);
nen = 3;
d = size(Dof,2);
dd = (d*d-d)/2;
X = reshape(Coord',nrdof,1);
Fext = zeros(nrdof,1);
FR = zeros(n,1);
FP = zeros(n,1);
ux_hist = zeros(n,1);
eps_p_bar = zeros(n,1);
x = X;
old_x = x;
bc0 = [dof_fixed,X(dof_fixed)];
old_state = struct('Fp',eye(d),'eps_p_bar',num2cell(zeros(nrelem,1)));
state = old_state;

ux = 0;
old_ux = 0;
old_dux = 1;
old_du = zeros(nrdof,1); % Used for initial guess for the increment (scaled increment)
i = 1;
conv_count = 0;
unload = false;

% Main loop
tic
while true
%%%%%%%%%%%%%%%% Special loading case, break condition part b
%    if old_ux >= umax
%        break;
%    end

%%%%%%%%%%%%%%%% Special loading case, part c
    if old_ux >= umax - min_dux && ~unload
        dux = -dux0;
        unload = true;
    end

    if old_ux + dux > umax
        dux = umax - old_ux;
    end
    ux = old_ux + dux;
    fprintf(' \t\t\tTimestep %4d, ux = %e, dux = %e, conv_count = %d\n',i,ux,dux,conv_count)
    bcx = [bc0;dof_prescr,X(dof_prescr)+ux]; % Boundary conditions on the absolute value
    x = old_x + old_du*(dux/old_dux); % Previous increment scaled
    x(bcx(:,1)) = bcx(:,2);
    conv = true;
    for p = 1:maxiter
        apos = 0; II = zeros(nrelem*(nen*d)^2,1); JJ = zeros(nrelem*(nen*d)^2,1); XX = zeros(nrelem*(nen*d)^2,1);
        T = zeros(nrdof,1);
        ex = extract(Edof,x);
        for j = 1:nrelem
            X1 = [Ex(j,:);Ey(j,:)]';
            x1 = reshape(ex(j,:),d,nen)';
            [d0N,dN] = cts(X1,x1);
            F = x1'*d0N';
%            [S,dS_dE] = neo_hooke(F,modpar,false);
%%{
            C = F'*F;
            %[S,state(j),econv] = neo_hooke_plast(C,old_state(j),modpar);
            [S,state(j),econv] = neo_hooke_visco(C,Dt,old_state(j),modpar);
            conv = econv && conv;
            % Numerical derivative
            h = 1e10*eps(max(C(:)));
            dS_dE = zeros(d*d);
            for k = 1:d*d
                Cd = m2v(C); Cd(k) = Cd(k)+h;
                %Sh = neo_hooke_plast(v2m(Cd),old_state(j),modpar);
                Sh = neo_hooke_visco(v2m(Cd),Dt,old_state(j),modpar);
                dS_dE(:,k) = 2/h*m2v(Sh-S);
            end
            dS_dE = [eye(d),zeros(d);zeros(d),ones(d)/2]*dS_dE;
%%}
            J = det(F);
            s = 1/J*F*S*F';
            c = 1/J*over_outer(F,F)*dS_dE*over_outer(F',F');

            Te = compute_T(x1,t,dN,s);
            [Kc,Ks] = compute_K(dN,x1,t,c,s);
            Ke = Kc+Ks;
            for ii = 1:nen*d, for jj = 1:nen*d
                apos = apos+1; 
                II(apos) = Edof(j,1+ii); JJ(apos) = Edof(j,1+jj); XX(apos) = Ke(ii,jj);
            end, end
            T(Edof(j,2:end)) = T(Edof(j,2:end)) + Te;
        end
        K = sparse(II,JJ,XX);

        bcu = [bcx(:,1),bcx(:,2)-x(bcx(:,1))];
        R = T-Fext;
        res = norm(R(dof_free));
        fprintf('Iteration %2d, Res = %e\n',p,res/norm(R(dof_fixed)))
        if res <= tol*norm(R(dof_fixed));
            break;
        end
        ddu = solveq(K,-R,bcu);
        x = x + ddu;
    end
    if p == maxiter || ~conv
        conv_count = 0;
        if abs(dux) < min_dux % then give up
            warning('No FE convergence')
            i = i-1;
            break;
        end
        dux = dux/2;
        continue;
    end
    conv_count = conv_count + 1;
    old_state = state;
    old_du = x - old_x;
    old_x = x;
    old_dux = dux;
    old_ux = ux;

    i = i + 1;
    FR(i) = sum(T(dof_vert));
    ux_hist(i) = ux;
    if conv_count > 2 && p <= maxiter/2 && ~unload
        dux = dux*1.5;
        if dux > max_dux
            dux = max_dux;
        end
    end

%%%%%%%%%%%%%%%% Special loading case, break condition for part c
    if unload && FR(i) >= 0
        break;
    end
end
toc

figure
plot(ux_hist(1:i)*1e3,-FR(1:i))
%fplot(@(e) interp1(ux_hist(1:i)*1e3,-FR(1:i),e),ux_hist([1,i])*1e3); % Somewhat optimized plot
ylabel('F_x\ \text{[N]}')
xlabel('u_x\ \text{[mm]}')
axis tight
matlab2tikz('ux_Fx2.tikz','width','\figurewidth','height','\figureheight')

%figure
%%plot(ux(1:i),eps_p_bar(1:i));

%figure
%plot(ux(1:i),Fp(1:i));

figure
%eldraw2(Ex,Ey,[1,1,0])
ed = extract(Edof,x-X);
eldisp2(Ex,Ey,ed,[1,2,0],1);
axis tight

figure
patch(ex(:,1:2:end)',ex(:,2:2:end)',[old_state.eps_p_bar]);
colorbar
colormap('hot')
axis equal
axis tight

