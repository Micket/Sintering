function [eps,sig,state,E_tang,lconv] = umaterial(oldeps,oldsig,state,Deps,Dt,modpar)
% [eps(:,i),sig(:,i),state(i),E_tang,lconv] = ...
%               umaterial(oldeps(:,i),oldsig(:,i),state(i),Deps,Dt,modpar);

E      = modpar.E;
H      = modpar.H;
nu     = modpar.nu;
sy     = modpar.sy;
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
eps_e_dev_tr = state.eps_e_dev + Idev*Deps;
eps_e_e_tr = sqrt(2/3*(eps_e_dev_tr'*eps_e_dev_tr));
sigma_e_tr = 3*G*eps_e_e_tr;
Phi_tr = sigma_e_tr - sy;
if Phi_tr <= 0 
    state.eps_e_dev = eps_e_dev_tr;
    sig = Ee*(eps_e_dev_tr+(eps'*I2)*I2)
    E_tang = Ee;
    lconv = 1;
else
    if n == 1
        eps_e_e = (sy*Dt + eps_e_e_tr*sc*t_star)/(Dt+sc*t_star);
        Phi = 3*G*eps_e_e - sy;
        bingham = Phi/sc;
        mu = Dt/t_star*bingham;
    else
        TOL = 1e-6*norm(Deps);
        itermax = 10;
        lconv = 0;
        eps_e_e = eps_e_e_tr;
        Jtmp = 3*G*n/sc;
        for i = 1:itermax
            Phi = 3*G*eps_e_e - sy;
            bingham = Phi/sc;
            mu = Dt/t_star*(bingham)^n;
            R = eps_e_e - eps_e_e_tr + mu;
            J = 1 + Jtmp*mu/bingham;
            eps_e_e = eps_e_e - R/J;
            if R < TOL || R == 0
                lconv = 1;
                break
            end
        end
        if lconv == 0
            warning(['Not converged:',num2str(R)]);
        end
    end
    eps_e_dev = eps_e_e/(eps_e_e+mu)*eps_e_dev_tr;
    state.eps_e_dev = eps_e_dev;
    eps_e = eps_e_dev + (eps'*I2)/3*I2;
    sig = Ee*eps_e;
    Q = Idev - 2/(3*eps_e_e^2)*eps_e_dev*eps_e_dev';
    ha = 3*G + t_star*sc/Dt*(sc/Phi)^(n-1);
    b = 1/(eps_e_e/mu + 1);
    E_tang = Ee - (2*G*b)*Q + -4*G^2/(ha*eps_e_e^2)*eps_e_dev*eps_e_dev';
end

