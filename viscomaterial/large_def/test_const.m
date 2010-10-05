clear all
close all

free = 0; % Run free deformation
%free = 1;
print = false;
maxiter = 1000; % Iterations for zero stress
nrtstep = 1000;
tmax = 9000;
tm = 3600; % Heat up time
thetamax = 1250;
theta0 = 20;

theta = @(t) interp1(...
    [0,tm,tmax],...
    [theta0,thetamax,thetamax],t); % Temperature ramp

d  = 2; % Number of dimensions to run
% Some constants
I = eye(d);
Iv = zeros(d*d,1); Iv(1:d) = 1;
Idev = (eye(d*d) - Iv*Iv'/d);
offdiag = (d*d-d)/2;

% Material parameters
rhop0 = 0.5;
E = 6e11;
nu = 0.3;
G = E/(2+2*nu);
Kb = E/(3-6*nu);
Lsolid = 1e-6;
gamma = 1;

mp.G = G;
mp.Kb = Kb;
mp.alpha = 0; %7e-6;
mp.B = 2/3*(3/4/pi)^(-1/3)*gamma/Lsolid;

mp.c = 0.1;
mp.n = 1;

mp.P = 1e6;
mp.tau_c = 1e8;

mp.tau_yd = 1.5e9;
mp.gammaT = 4.6e4;
mp.thetaT = 601;

mp.t_stard = 8.6e-3;
mp.gammat = 2.0e4;
mp.thetat = 993;

% Not quite material parameter, but it doesn't matter.
mp.theta0 = 20;

state = struct('iFp',eye(d),'rhop',rhop0,'rhop0',rhop0);

% Basic check, how long time should this roughly take to go from rhop = 0.5 -> 0.6 ?
%t_star_approx = mp.t_stard*exp(mp.gammat/thetamax*(1-thetamax/mp.thetat));
%tau_y_approx = mp.tau_yd*exp(mp.gammaT/thetamax*(1-thetamax/mp.thetaT));
%a2_approx = mp.c2*(1-0.5)^mp.n2;
%Mmbh_approx = 3e6;
%Phi_approx =  1/mp.P * (a2_approx*Mmbh_approx^2 - tau_y_approx^2);
%eta_approx = mp.eta0*(Phi_approx/mp.tau_c)^mp.n;
%Dt_approx = log(0.6/0.5)*t_star_approx*3*mp.G/(2*a2_approx*eta_approx*Mmbh_approx)
%pause(2)

% Values to store for each iteration
sigma = zeros(d*d,nrtstep);
%k = zeros(1,nrtstep);
rhop = zeros(1,nrtstep); rhop(1) = rhop0;
plast = zeros(1,nrtstep); plast(1) = 1;
strain = zeros(1,nrtstep);

t = linspace(0,tmax,nrtstep);
if ~free 
    % Fixed 0 deformation
    C_old = eye(d);
    for i=2:nrtstep
        Dt = t(i)-t(i-1);
        F = I;
        C = F'*F;
        [S2,state,conv] = sintering_full(C_old,C,Dt,theta(t(i)),state,mp,true);
        C_old = C;
        EE = 1/2*(C-I);
        J = det(F);
        sigma(:,i) = m2v(1/J*F*S2*F');
        %k(i) = state.k;
        rhop(i) = state.rhop;
        plast(i) = 1/det(state.iFp);
        strain(i) = EE(1,1);
        if ~conv 
            break
            i = i-1;
        end
    end
else
    % Free sintering
    C = eye(d);
    C_new = C;
    for i=2:nrtstep
        fprintf(' Timestep %4d/%d, t = %e, theta = %4.0f, rhop = %e, C_11 = %e\n',...
            i,nrtstep,t(i),theta(t(i)),state.rhop,C(1,1));
        Dt = t(i)-t(i-1);
        for p = 1:maxiter
            [S,state_new,La,conv] = sintering_full(C,C_new,Dt,theta(t(i)),state,mp,print);
            if ~conv
                return
            end
            R = m2v(S-eye(d)*0e6);
            fprintf('Iteration %2d, ',p)
            fprintf('Stress = '), fprintf('%13.6e, ',S), fprintf('| Strain = '), fprintf('%13.6e, ',C), fprintf('\n')
            conv_stress = norm(R) < 10;
            if conv_stress || ~conv
                break
            end

            %dCv = -La\R;
            dCv = -La(1:d+offdiag,1:d+offdiag)\R(1:d+offdiag);
            dCv = [dCv;dCv(1+d:end)];
            C_new = C_new + v2m(dCv);
        end
        if ~conv
            warning('Material routine did not converge')
            i = i -1; break;
        end
        if ~conv_stress
            warning('Stress not converged')
            i = i -1; break;
        end

        C = C_new;
        EE = 1/2*(C-I);
        J = sqrt(det(C));
        % Stored values
        %sigma(:,i) = m2v(1/J*F*S*F');
        sigma(:,i) = m2v(S);
        %k(i) = state.k;
        rhop(i) = state.rhop;
        plast(i) = 1/det(state.iFp);
        strain(i) = EE(1,1);
        state = state_new;
        if ~conv 
            break
            i = i-1;
        end
    end
end

sm   = Iv'*sigma;
sdev = Idev*sigma;
se   = sqrt(3/2*sum(sdev.^2));

figure
fplot(@(e) interp1(t(1:i),rhop(1:i),e),t([1,end]),'r')
ylim([0,1])
xlabel('t\ \text{[s]}')
ylabel('\rho''\ \text{[-]}')
title('Relative density')


%figure
%fplot(theta,t([1,end]))
%xlabel('t\ \text{[s]}')
%ylabel('\theta\ \text{[C$^\circ$]}')
%title('Temp')

if ~free
figure, hold on
fplot(@(e) interp1(t(1:i),sigma(1,1:i),e),t([1,end]),'k--')
fplot(@(e) interp1(t(1:i),sigma(2,1:i),e),t([1,end]),'k-.')
fplot(@(e) interp1(t(1:i),sigma(3,1:i),e),t([1,end]),'k:')
fplot(@(e) interp1(t(1:i),se(1:i),e),t([1,end]),'b-')
fplot(@(e) interp1(t(1:i),sm(1:i),e),t([1,end]),'r-')
plot(t([1,end]),mp.tau_yd*[1,1],'g-')
xlabel('t\ \text{[s]}')
ylabel('\sigma\ \text{[Pa]}')
legend('\sigma_{11}','\sigma_{22}','\sigma_{12}','\sigma_e','\sigma_m')
%matlab2tikz('strain_stress.tikz','height','\figureheight','width','\figurewidth')
%matfig2pgf('filename','strain_stress.pgf','figwidth','\figurewidth')
end

figure
fplot(@(e) interp1(t(1:i),plast(1:i),e),t([1,end]),'r')
xlabel('t\ \text{[s]}')
title('Plastic deformation')

figure
fplot(@(e) interp1(t(1:i),strain(1:i),e),t([1,end]),'r')
xlabel('t\ \text{[s]}')
ylabel('E_{11}\ \text{[-]}')
title('Strain')

