function [eps,sig,state,E_tang,lconv] = umaterial(oldeps,oldsig,state,Deps,Dt,modpar)
% [eps(:,i),sig(:,i),state(i),E_tang,lconv] = ...
%               umaterial(oldeps(:,i),oldsig(:,i),state(i),Deps,Dt,modpar);

E      = modpar.E;
H      = modpar.H;
nu     = modpar.nu;
sc     = modpar.sc;
t_star = modpar.t_star;
n      = modpar.n;

L = E*nu/((1+nu)*(1-2*nu));
G = E/(2+2*nu);

I2 = zeros(9,1); I2(1:3) = 1; 
I = eye(9);
Idev = I - 1/3*I2*I2';
Ee = 2*G*I + L*I2*I2';

eps = oldeps + Deps;
eps_e_tr = state.eps_e + Deps;
sigma_tr = Ee*eps_e_tr;
sigma_dev_tr = Idev*sigma_tr;
sigma_e_tr = sqrt(3/2*(sigma_dev_tr'*sigma_dev_tr));
Phi_tr = sigma_e_tr;
if Phi_tr == 0 
    sig = sigma_tr;
    state.eps_e = eps_e_tr;
    E_tang = Ee;
    lconv = 1;
else
    TOL = 1e-9;
    itermax = 10;
    lconv = 0;
    eps_e = eps_e_tr;
    for i = 1:itermax
        sigma = Ee*eps_e;
        sigma_dev = Idev*sigma;
        sigma_e = sqrt(3/2*(sigma_dev'*sigma_dev));
        Phi = sigma_e;
        mu = (Phi/sc)^n*Dt/t_star;

        % Assemble R        
        Re = eps_e - eps_e_tr + mu*3/(2*sigma_e)*sigma_dev;
        Ru = (Phi/sc)^n- t_star/Dt*mu;
        R = [Re;Ru];
       
        % Assemble J
        Jee = I + 3*mu*G/sigma_e*(Idev-(3/(2*sigma_e^2)*sigma_dev)*(sigma_dev'));
        Jeu = 3/(2*sigma_e)*sigma_dev;
        Jue = 3*G*n/(sc*sigma_e)*(Phi/sc)^(n-1)*sigma_dev';
        Juu = -t_star/Dt;
        J = [Jee,Jeu;...
             Jue,Juu];

        if norm(Re) <= TOL*norm(Deps) && abs(Ru) < TOL
            lconv = 1;
            break
        end
        dX = J\R;
        eps_e = eps_e - dX(1:9);
        mu    = mu    - dX(10);
    end
    if lconv == 0
        (eps_e - eps_e_tr)/norm(Deps)
        [norm(Re)/norm(Deps),abs(mu)]
    end
    state.eps_e = eps_e;
    sig = sigma;
    Ji = inv(J);
    epart = 1:9;
    E_tang = Ee*Ji(epart,epart);
end

