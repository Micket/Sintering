function [eps,sig,state,E_tang,lconv] = umaterial(oldeps,oldsig,state,Deps,Dt,modpar)
% [eps(:,i),sig(:,i),state(i),E_tang,lconv] = ...
%               umaterial(oldeps(:,i),oldsig(:,i),state(i),Deps,Dt,modpar);

E      = modpar.E;
H      = modpar.H;
nu     = modpar.nu;
sy     = modpar.sy;
sy_dyn = modpar.sy_dyn;
t_star = modpar.t_star;
kappa_inf = modpar.kappa_inf;
n      = modpar.n;
S      = modpar.S;
m      = modpar.m;
b      = state.b;
k      = state.k;

L = E*nu/((1+nu)*(1-2*nu));
G = E/(2+2*nu);

I2 = zeros(9,1); I2(1:3) = 1; 
I = eye(9);
Idev = I - 1/3*I2*I2';
Ee = 2*G*I + L*I2*I2';

eps = oldeps + Deps;
eps_e_tr = state.eps_e + Deps;
sigma = Ee*eps_e_tr;
sigma_dev = Idev*sigma;
sigma_e = sqrt(3/2*(sigma_dev'*sigma_dev));
kappa = -H*k;
if sigma_e - (sy + kappa) <= 0
    sig = sigma;
    state.eps_e = eps_e_tr;
    E_tang = Ee;
    lconv = 1;
else
    TOL = 1e-9;
    c = sy_dyn - sy;
    itermax = 20;
    lconv = 0;
    oldk = state.k;
    mu = state.mu;
    k = oldk;
    eps_e = eps_e_tr;
    for i = 1:itermax
        kappa = -H*k;
        sigma = Ee*eps_e;
        sigma_dev = Idev*sigma;
        sigma_e = sqrt(3/2*(sigma_dev'*sigma_dev));
        Phi = sigma_e - (sy + kappa);
        beta = -(eps_e'*sigma+H*k^2)/2;
        q = Phi/(c-Phi);
        
        % Assemble J
        Jee = I + 3*mu*G/sigma_e*(Idev-(3/(2*sigma_e^2)*sigma_dev)*(sigma_dev'));
        Jek = zeros(9,1);
        Jeu = 3/(2*sigma_e)*sigma_dev;
        
        Jke = zeros(1,9);
        Jkk = 1+mu*H/kappa_inf; 
        Jku = 1-kappa/kappa_inf;
    
        if q > 0
            Jtmp = n*q^(n-1)*c/(c-Phi)^2;
            Jue = (Jtmp*3*G/sigma_e)*sigma_dev';
            Juk = Jtmp*H;
        else
            Jue = zeros(1,9);
            Juk = 0;
        end
        Juu = -t_star/Dt;

        J = [Jee,Jek,Jeu;...
             Jke,Jkk,Jku;...
             Jue,Juk,Juu];

        % Assemble R        
        Re = eps_e-eps_e_tr+mu*3/(2*sigma_e)*sigma_dev;
        Rk = k - oldk + mu*(1-kappa/kappa_inf);
        Ru = -t_star/Dt*mu;
        if q > 0
            Ru = Ru + q^n;
        end
        R = [Re;Rk;Ru];

        if norm(R) < TOL
            lconv = 1;
            break
        end

        dX = -J\R;

        eps_e = eps_e + dX(1:9);
        k     = k     + dX(10);
        mu  = mu      + dX(11);
    end
    state.eps_e = eps_e;
    state.k = k;
    state.mu = mu;
    sig = sigma;
    Ee = Ee;
    Ji = inv(J);
    epart = 1:9;
    E_tang = Ee*Ji(epart,epart);
end

