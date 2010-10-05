function [S2,state,conv] = neo_hook_plast(C,state,modpar)
mu = modpar.mu;
lambda = modpar.lambda;
H = modpar.H;
tau_y = modpar.tau_y;

maxiter = 10;
tol = 1e-9;
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
Phi = sqrt(3/2)*norm(Mdev(:))-(tau_y + H*state.eps_p_bar);
if Phi <= 0
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
        R = sqrt(3/2)*norm(Mdev(:))-(tau_y + H*(state.eps_p_bar+Dy));
        if abs(R) < tol*tau_y
            break
        end
        %h = 1e-10;
        %J = (RR(iFp_old,nu_old,C,eps_p_bar,modpar,h) - RR(iFp_old,nu_old,C,eps_p_bar,modpar,0))/h;

        % Analytical derivative
        dMdCe = over_outer(I,S2t)+lambda/2*outer(I,iCe)+(mu-lambda*log(Je))*over_outer(I,iCe);
        dCedDy = m2v(-2*iFp'*C*iFp_old*nu_old);
        nu_v = m2v(sqrt(3/2)*Mdev/norm(Mdev(:)));
        J = nu_v'*dMdCe*dCedDy-H;
        Dy = Dy - R/J;
    end
    conv = abs(R) < tol*tau_y;
    if ~conv
        warning(['Error not converged, res = ',num2str(abs(R)/tau_y)])
    end
    S2 = iFp*S2t*iFp';
    state.Fp = inv(iFp);
    state.eps_p_bar = eps_p_bar + Dy;
end
end

%{
function R = RR(iFp_old,nu_old,C,eps_p_bar,modpar,Dy)
    d = size(C,1);
    I = eye(d);
    iFp = iFp_old - Dy*iFp_old*nu_old;
    Ce = iFp'*C*iFp;
    Je = sqrt(det(Ce));
    iCe = inv(Ce);
    S2t = modpar.mu*(I-iCe)+modpar.lambda*log(Je)*iCe;
    M = Ce*S2t;
    Mdev = M - trace(M)/d*I;
    R = sqrt(3/2)*norm(Mdev(:))-(modpar.tau_y + modpar.H*(eps_p_bar+Dy));
end
%}