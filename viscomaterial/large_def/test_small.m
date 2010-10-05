clear all
%close all

free = 0; % Run free deformation
free = 1;
print = false;
%print = true;
maxiter = 20; % Iterations for zero stress
nrtstep = 1000;
tmax = 900000;
tm = 360000; % Heat up time
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
rhop0 = 0.6;
%E = 6e11;
E = 1e9;
nu = 0.3;
G = E/(2+2*nu);
Kb = E/(3-6*nu);
Lsolid = 1e-6;
gamma = 1;

mp.G = G;
mp.Kb = Kb;
mp.alpha = 7e-6;
mp.B = 2/3*(3/4/pi)^(-1/3)*gamma/Lsolid;

mp.c = 0.1;
mp.n = 1;

mp.P = 1e7; % Just leaving this at 1 MPa, not pretty
mp.S_c = 5e11;

mp.S_yd = 1.5e9;
mp.gammaT = 4.6e4;
mp.thetaT = 601;

mp.t_stard = 8.6e-3;
mp.gammat = 2.0e4;
mp.thetat = 993;

% Not quite material parameter, but it doesn't matter.
mp.theta0 = 20;

state = struct('E_p',zeros(d),'rhop',rhop0,'rhop0',rhop0);

% Basic check, how long time should this roughly take to go from rhop = 0.5 -> 0.6 ?
% This doens't apply anymore
%t_star_approx = mp.t_stard*exp(mp.gammat/thetamax*(1-thetamax/mp.thetat));
%S_y_approx = mp.S_yd*exp(mp.gammaT/thetamax*(1-thetamax/mp.thetaT));
%a_approx = mp.c*(1-0.5);
%Smh_approx = 3e6;
%Phi_approx =  1/mp.P * (a_approx*Smh_approx^2 - S_y_approx^2);
%eta_approx = (Phi_approx/mp.S_c)^mp.n;
%Dt_approx = log(0.6/0.5)*t_star_approx*mp.P/(2*a_approx*eta_approx*Smh_approx)
%pause(2)

% Values to store for each iteration
sigma = zeros(d*d,nrtstep);
rhop = zeros(1,nrtstep); rhop(1) = rhop0;
plast = zeros(1,nrtstep);
strain = zeros(1,nrtstep);

t = linspace(0,tmax,nrtstep);
if ~free 
    % Fixed 0 deformation
    E_old = zeros(d);
    for i=2:nrtstep
        Dt = t(i)-t(i-1);
        E = zeros(d);
        [S,state,conv] = sintering_small(E_old,E,Dt,theta(t(i-1)),theta(t(i)),state,mp,true);
        sigma(:,i) = m2v(S);
        rhop(i) = state.rhop;
        plast(i) = trace(state.E_p);
        strain(i) = E(1,1);
        if ~conv 
            break
            i = i-1;
        end
    end
else
    % Free sintering
    E = zeros(d);
    E_new = E;
    tic
    for i=2:nrtstep
        fprintf(' Timestep %4d/%d, t = %e, theta = %4.0f, rhop = %e, E_11 = %e\n',...
            i,nrtstep,t(i),theta(t(i)),state.rhop,E(1,1));
        Dt = t(i)-t(i-1);
        Dtheta = theta(t(i)) - theta(t(i-1));
        DE_free = I*(mp.alpha*Dtheta);
        %E_new = E_new + DE_free;
        for p = 1:maxiter
            [S,state_new,Ea,conv] = sintering_small(E,E_new,Dt,theta(t(i-1)),theta(t(i)),state,mp,print);
            R = m2v(S); % Zero stress!
            fprintf('Iteration %2d, ',p)
            fprintf('Stress = '), fprintf('%13.6e, ',S), fprintf('| Strain = '), fprintf('%13.6e, ',E_new), fprintf('\n')
            conv_stress = norm(R) < 10;
            if conv_stress || ~conv
                break
            end

            dEv = -Ea\R;
            E_new = E_new + v2m(dEv);
        end
        if ~conv
            warning('Material routine did not converge')
            i = i -1;
            return
        end
        if ~conv_stress
            warning('Stress not converged')
            i = i -1;
            return;
        end

        E = E_new;
        sigma(:,i) = m2v(S);
        rhop(i) = state.rhop;
        plast(i) = trace(state.E_p)/d;
        strain(i) = E(1,1);
        state = state_new;
        if ~conv 
            break
            i = i-1;
        end
    end
    toc
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

if ~free
    figure, hold on
    fplot(@(e) interp1(t(1:i),sigma(1,1:i),e),t([1,end]),'k--')
    fplot(@(e) interp1(t(1:i),sigma(2,1:i),e),t([1,end]),'k-.')
    fplot(@(e) interp1(t(1:i),sigma(3,1:i),e),t([1,end]),'k:')
    fplot(@(e) interp1(t(1:i),se(1:i),e),t([1,end]),'b-')
    fplot(@(e) interp1(t(1:i),sm(1:i),e),t([1,end]),'r-')
    plot(t([1,end]),mp.S_yd*[1,1],'g-')
    xlabel('t\ \text{[s]}')
    ylabel('\sigma\ \text{[Pa]}')
    legend('\sigma_{11}','\sigma_{22}','\sigma_{12}','\sigma_e','\sigma_m')
    %matlab2tikz('strain_stress.tikz','height','\figureheight','width','\figurewidth')
    %matfig2pgf('filename','strain_stress.pgf','figwidth','\figurewidth')
else

end

figure, hold on
fplot(@(e) interp1(t(1:i),strain(1:i),e),t([1,end]),'r')
fplot(@(e) interp1(t(1:i), plast(1:i),e),t([1,end]),'b')
%ylim([-0.15*2,0.05])
legend('Real','Plastic')
xlabel('t\ \text{[s]}')
title('Strain')

