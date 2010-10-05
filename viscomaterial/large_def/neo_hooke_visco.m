function [S2,state,conv] = neo_hook_visco(C,Dt,state,modpar)
%function [S2,state,conv] = neo_hook_visco(C,Dt,state,modpar)

mu     = modpar.mu;
lambda = modpar.lambda;
H      = modpar.H;
tau_y  = modpar.tau_y;
tau_c  = modpar.tau_c;
t_star = modpar.t_star;
eta0   = modpar.eta0;
n      = modpar.n;
c      = modpar.c;

maxiter = 10;
tol = 1e-16;
d = size(C,1);
I = eye(d);

Fp = state.Fp;
eps_p_bar = state.eps_p_bar;

Ce = Fp'\C/Fp;
Je = sqrt(det(Ce));
iCe = inv(Ce);
S2t = mu*(I-iCe)+lambda*log(Je)*iCe;
M = Ce*S2t;
Mdev = M - trace(M)/d*I;
Phi = sqrt(3/2)*norm(Mdev(:))-(tau_y + H*eps_p_bar);
if Phi <= 0
    disp('Elast')
    S2 = Fp\S2t/Fp';
    conv = true;
else
    nu_old = sqrt(3/2)*Mdev/norm(Mdev(:));
    nu_old_v = m2v(nu_old);
    Fp_old = Fp;
    iFp_old = inv(Fp_old);
    Dy = 0;
    for k = 1:maxiter
        iFp = iFp_old - Dy*iFp_old*nu_old;
        Ce = iFp'*C*iFp;
        Je = sqrt(det(Ce));
        iCe = inv(Ce);
        S2t = mu*(I-iCe)+lambda*log(Je)*iCe;
        M = Ce*S2t;
        Mdev = M - trace(M)/d*I;
        Phi = sqrt(3/2)*norm(Mdev(:))-(tau_y + H*(eps_p_bar+Dy));
        eta = 0;
        if Phi > 0
            %eta = eta0*(Phi/(c-Phi))^n;
            eta = eta0*(Phi/tau_c)^n;
        end
        R = eta*Dt - t_star*Dy;
        %R = RR(Dt,iFp_old,nu_old,C,eps_p_bar,modpar,Dy);
        %fprintf('Res %d, %e\n',k,abs(R));
        if abs(R) < tol*mu
            break
        end
        %h = 1e-10;
        %J1 = (RR(Dt,iFp_old,nu_old,C,eps_p_bar,modpar,Dy+h) - R)/h;

        % Analytical derivative
        detadDy = 0;
        if Phi > 0 % if iterations might to back into a elastic region.
            dMdCe = over_outer(I,S2t)+lambda/2*outer(I,iCe)+(mu-lambda*log(Je))*over_outer(I,iCe);
            dCedDy = m2v(-2*iFp'*C*iFp_old*nu_old);
            nu_v = m2v(sqrt(3/2)*Mdev/norm(Mdev(:)));
            dPhidDy = nu_v'*dMdCe*dCedDy-H;
            %detadDy = eta0*n*(Phi/(c-Phi))^(n-1)*c/(c-Phi)^2*dPhidDy;
            detadDy = eta0*n*(Phi/tau_c)^(n-1)/tau_c*dPhidDy;
        end
        J = detadDy*Dt - t_star;
        %dJ = abs((J-J1)/J)

        Dy = Dy - R/J;
    end
    fprintf('iterations = %2d, Dy = %e\n',k,Dy)
    conv = abs(R) < tol*mu;
    if ~conv
        warning(['Error not converged, res = ',num2str(abs(R)/mu)])
    end
    iFp = iFp_old - Dy*iFp_old*nu_old;
    S2 = iFp*S2t*iFp';
    state.Fp = inv(iFp);
    state.eps_p_bar = eps_p_bar + Dy;
end
end

%{
function R = RR(Dt,iFp_old,nu_old,C,eps_p_bar,modpar,Dy)
    d = size(C,1);
    I = eye(d);
    iFp = iFp_old - Dy*iFp_old*nu_old;
    Ce = iFp'*C*iFp;
    Je = sqrt(det(Ce));
    iCe = inv(Ce);
    S2t = modpar.mu*(I-iCe)+modpar.lambda*log(Je)*iCe;
    M = Ce*S2t;
    Mdev = M - trace(M)/d*I;
    Phi = sqrt(3/2)*norm(Mdev(:))-(modpar.tau_y + modpar.H*(eps_p_bar+Dy));
    eta = 0;
    if Phi > 0
        %eta = modpar.eta0*(Phi/(modpar.c-Phi))^modpar.n;
        eta = modpar.eta0*(Phi/modpar.tau_c)^modpar.n;
    end
    R = eta*Dt - modpar.t_star*Dy;
end
%}
