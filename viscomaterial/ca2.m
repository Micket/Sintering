clear all
%close all

col = 'rbkmc';

% Material data
modpar.E = 150e9;
modpar.H = 150e9;
modpar.nu = 0.3;
modpar.sy = 0e6;
modpar.sc = 10e6;
modpar.n = 4;
modpar.t_star = 1e4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('a) Performance of the driver')
Dt = [1e-1,1e-0,1e+1,1e+2]/1000;
Deps_p = cell(length(Dt),1);
for j = 1:length(Dt)
    Deps = [1e-5,0,0,0,0,0,0,0,0]'
    steps = round(0.001/Deps(1));
    eps = zeros(9,steps); sig = zeros(9,steps);
    state = struct('eps_e',zeros(9,1),'eps_e_dev',zeros(9,1),'mu_h',num2cell(zeros(steps,1)));

    for i = 2:steps
        [eps(:,i),sig(:,i),state(i),E_tang,lconv] = viscored(eps(:,i-1),sig(:,i-1),state(i-1),Deps,Dt(j),modpar);
        if lconv < 1
            fprintf('Not converged\n')
        end
    end
    Idev = eye(9); Idev = Idev-[1,1,1,0,0,0,0,0,0]'*[1,1,1,0,0,0,0,0,0]/3;
    sig_dev = Idev*sig; 
    sig_e = sqrt(3/2*sum(sig_dev.^2,1));

    subplot(2,1,1), plot(eps(1,1:i),sig_e(1:i),col(j)), hold on
    xlabel('$\epsilon_{11}$','interpreter','latex')
    ylabel('$\sigma_e$','interpreter','latex')
    subplot(2,1,2), plot(eps(1,1:i),sig(1,1:i),col(j)), hold on
    xlabel('$\epsilon_{11}$','interpreter','latex')
    ylabel('$\sigma_{11}$','interpreter','latex')
    Deps_p{j} = sprintf('$\\dot{\\epsilon}_{11}=%.2e$',Deps(1)/Dt(j));
end
figure(1), H=legend(Deps_p), set(H,'interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('b) Creep behaviour')
sigma_0 = 2*modpar.sc;
sigma_11 = @(t) sigma_0*t/50+sigma_0*(1-t/50)*heaviside(t-50+1e-5);
maxsteps = 10000;
eps = zeros(9,maxsteps); sig = zeros(9,maxsteps);
state = struct('eps_e',zeros(9,1),'eps_e_dev',zeros(9,1),'mu_h',num2cell(zeros(maxsteps,1)));
time = zeros(1,maxsteps);
TOL = 1; itermax = 10;
Dt = 0.01;
Dtmin = 0.0001;
Deps = [0,0,0,0,0,0,0,0,0]';
for i = 2:maxsteps
    time(i) = time(i-1) + Dt;
    sigma = [sigma_11(time(i));0;0;0;0;0;0;0;0];
    for j=1:itermax
        [eps(:,i),sig(:,i),state(i),E_tang,lconv] = viscored(eps(:,i-1),sig(:,i-1),state(i-1),Deps,Dt,modpar);
        rB = sigma - sig(:,i);
        if max(abs(rB)) < TOL
            break
        end
        Deps = Deps + E_tang\rB;
    end
    if j == itermax
        error('Stress not converged\n',i)
    end
    fprintf('time = %f, eps_11^p = %1.2e\n',time(i),Deps(1)/Dt);
end

figure(2)
subplot(2,1,1),fplot(sigma_11,[0,time(i)])
xlabel('$t$ [t]','interpreter','latex','fontsize',18)
ylabel('$\sigma_{11}$ [Pa]','interpreter','latex','fontsize',18)
subplot(2,1,2),plot(time(1:i),eps(1,1:i)), hold on
xlabel('$t$ [s]','interpreter','latex','fontsize',18)
ylabel('$\epsilon$','interpreter','latex','fontsize',18)

