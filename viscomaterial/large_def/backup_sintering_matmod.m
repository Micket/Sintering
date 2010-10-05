function [S2,state,conv] = sintering_matmod(C_old,C,Dt,theta,state,mp,print)
%function [S2,state,conv] = sintering_matmod(C_old,C,Dt,theta,state,mp,print)
d = size(C,1);
I = eye(d);

maxiter = 50;
E = 1/2*(C-I);
tolerance = min([1e-12,norm(E)*1e-6]);

tau_y  = mp.tau_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
%H      = mp.Hd     *exptemp_fun(mp.gammaH,theta,mp.thetaH);
f      = 1+mp.alpha*(theta-mp.theta0);

Fp_old   = state.Fp;
k_old    = state.k;
rhop_old = state.rhop;
rhop0    = state.rhop0;

Cb_old = (Fp_old')\C_old/Fp_old;
%eb_old = logm(Cb_old)/2;
eb_old = 1/2*(Cb_old - I);
eb_old_vol = trace(eb_old);
eb_old_dev = eb_old - eb_old_vol/d*I;
eb_old_e = sqrt(2/3)*norm(eb_old_dev);

Cb_tr = (Fp_old')\C/Fp_old;
%eb_tr = logm(Cb_tr)/2;
eb_tr = 1/2*(Cb_tr - I);
eb_tr_vol = trace(eb_tr);
eb_tr_dev = eb_tr - eb_tr_vol/d*I;
eb_tr_e = sqrt(2/3)*norm(eb_tr_dev);

X = [rhop_old,0,eb_old_e,eb_old_vol]';

for j = 1:maxiter
    % Analytical jacobian
    %k = X(1); rhop = X(2); Dy = X(3); eb_e = X(4); eb_vol = X(5);
    rhop = X(1); Dy = X(2); eb_e = X(3); eb_vol = X(4);
    %ebb_vol = eb_vol - d*log(f);
    ebb_vol = eb_vol - d*log(f);
    
    %K = -H*k; 
    %Kd = -mp.Hd*k;
    % Common expressions (and derivatives w.r.t. k and rhop)
    %tau = (tau_y/(tau_y+K))^2; taup = (2*H*tau_y^2)/(tau_y + K)^3;
    %taud = (mp.tau_yd/(mp.tau_yd-mp.Hd*k))^2; taudp = (2*mp.Hd*mp.tau_yd^2)/(mp.tau_yd + Kd)^3;
    q = rhop0/((1-rhop)*rhop^2)^(1/3); qp  = rhop0*(3*rhop-2)/(3*(1-rhop)^(4/3)*rhop^(5/3));
    a1 = 1+(mp.c1-1)*(1-rhop)^mp.n1;   a1p = -(mp.c1-1)*mp.n1*(1-rhop)^(mp.n1-1);
    a2 = mp.c2*(1-rhop)^mp.n2;         a2p = -mp.c2*mp.n2*(1-rhop)^(mp.n2-1);

    Meb = eb_e*3*mp.G;
    Mmb = ebb_vol*mp.Kb;
    Mmbh = Mmb - mp.B*q;

    Phi  = 1/(3*mp.G) * (a1*Meb^2  + a2*Mmbh^2 - tau_y^2); % *tau
    %Phid = 1/(3*mp.G) * (a1*taud*Meb^2 + a2*Mmbh^2 - mp.tau_yd^2);
    Phid = -mp.tau_yd^2/(3*mp.G);
    if Phid > 0
        Meb
        Mmb
        Mmbh
        keyboard
        error('Inaccessible state! This shouldn''t happen.')
    end
    eta = 0;
    etap = 0;
    if Phi > 0 
        eta = mp.eta0*(-Phi/Phid)^mp.n;
        etap = mp.eta0*mp.n*(-Phi/Phid)^(mp.n-1);
    end

    % Derivatives of the yield surfaces.
    %dPdk  = 1/(3*mp.G) * Meb^2*a1*taup;
    %dPddk = 1/(3*mp.G) * Meb^2*a1*taudp;
    %dPddk = 0;
    %dPdrhop  = 1/(3*mp.G) * (a1p*tau*Meb^2   + (a2p*Mmbh^2 - 2*a2*Mmbh*mp.B*qp));
    dPdrhop  = 1/(3*mp.G) * (a1p*Meb^2   + (a2p*Mmbh^2 - 2*a2*Mmbh*mp.B*qp));
    %dPddrhop = 1/(3*mp.G) * (a1p*taud*Meb^2  + (a2p*Mmbh^2 - 2*a2*Mmbh*mp.B*qp));
    dPddrhop = 0;
    %dPdMeb = 2*a1*eb_e*tau;
    dPdMeb = 2*a1*eb_e;
    %dPddMeb = 2*a1*eb_e*taud;
    dPddMeb = 0;
    dPdMmbh = 2*a2*Mmbh/(3*mp.G);
    %dPddMmbh = dPdMmbh; 
    dPddMmbh = 0; 
    dPdMmb = dPdMmbh;
    %dPddMmb = dPdMmb;
    dPddMmb = 0;

    % Residuals
    %Rk = k - k_old + Dy*dPdMeb;
    Rp = log(rhop) - log(rhop_old) + Dy*dPdMmbh;
    Rg = t_star*Dy - Dt*eta;
    Re = eb_e - eb_tr_e + Dy*dPdMeb;
    Rm = eb_vol - eb_tr_vol + Dy*dPdMmbh;
    R = [Rp;Rg;Re;Rm];
    %R = [Rk;Rp;Rg;Re;Rm];

    if print
        fprintf('Constitutive iteration %2d, Res = ',j); fprintf('%11.4e, ',R); fprintf('\n');
        %fprintf('                           X   = ',j); fprintf('%11.4e, ',X); fprintf('\n');
    end
    conv = all(abs(R) < tolerance);
    if conv
        break
    end

    % Jacobian
    %Jk = [1 + 2*Dy*a1*eb_e*taup,...
    %      2*Dy*a1p*eb_e*tau,...
    %      dPdMeb,...
    %      2*Dy*a1*tau,...
    %      0];
    
    Jp = [...%0,...
          1/rhop + 2*Dy/(3*mp.G)*(a2p*Mmbh - a2*mp.B*qp),...
          dPdMmbh,...
          0,...
          Dy*2*a2*mp.Kb/(3*mp.G)];

    h = Dt*etap/Phid^2;
    Jg = [...%h*(Phid*dPdk - Phi*dPddk),...
          h*(Phid*dPdrhop - Phi*dPddrhop),...
          t_star,... % Done
          h*3*mp.G*(Phid*dPdMeb - Phi*dPddMeb),...
          h*mp.Kb*(Phid*dPdMmb - Phi*dPddMmb)];

    Je = [...%2*Dy*a1*eb_e*taup,...
          2*Dy*a1p*eb_e,... % *tau
          dPdMeb,...
          1+2*Dy*a1,... % *tau
          0];

    Jm = [...%0,...
          2*Dy/(3*mp.G)*(a2p*Mmbh-a2*mp.B*qp),...
          dPdMmbh,...
          0,...
          1+Dy*2*a2*mp.Kb/(3*mp.G)];
          
    J = [Jp; Jg; Je; Jm]; % Jk

    % Debugging; numerical jacobian
    %{
    [Rn,Pn,Pdn] = Res2(mp,theta,Dt,rhop0,k_old,rhop_old,eb_tr_e,eb_tr_vol,X);
    hn = 1e-4*[1e-6,1e-6,1e-6,1e-6,1e-6]';
    Jn = zeros(length(hn));
    for i = 1:length(hn);
        Xh = X; Xh(i) = Xh(i) + hn(i);
        Rh = Res2(mp,theta,Dt,rhop0,k_old,rhop_old,eb_tr_e,eb_tr_vol,Xh);
        Jn(:,i) = (Rh-Rn)/hn(i);
    end
    % Check for any difference
    dJ = Jn-J;
    dR = Rn-R;
    if norm(dJ(:))/norm(J(:)) > 1e-0 || norm(dR) > 1e-1
        dJ
        keyboard
    end
    %}

    %if rcond(J) < 1e-16
    %    keyboard
    %end
    %X = X -Jn\Rn;
    X = X -J\R;
end

if print
    fprintf('Constitutive iterations = %2d, Dy = %e\n',j,Dy)
end
if ~conv
    warning(['Material not converged, res = ',num2str(R')])
    conv = 1;
end
% Postprocessing, calculating the full forms
if eb_e == 0
    eb_dev = zeros(d)*eb_tr_dev;
else
    eb_dev = eb_e/(eb_e + Dy*dPdMeb)*eb_tr_dev;
end
Mb = 2*mp.G*eb_dev + Mmb*I;
nu = 2*a1*eb_dev + 1/d*dPdMmbh*I;
Fp = expm(Dy*nu)*Fp_old;
Cb = (Fp')\C/Fp;

%eb_vol2 = log(det(Cb))/d;
%(eb_vol-eb_vol2)/eb_vol

S2b = Mb; % it's actually Cb\Mb, but since i changed free energy, i work directly with S2 (but it's just a name)
S2 = Fp\S2b/(Fp');

state.Fp   = Fp;
state.rhop = rhop;
%state.k    = k;

% ATS Tensor (not done, dA^(-1) is missing.)
I4 = eye(d*d);
Iv = zeros(d*d,1); Iv(1:d) = 1;
% Elastic part
L2e = 2*mp.G*(I4-1/3*Iv'*Iv) + mp.Kb*Iv*Iv';
% Plastic part
%iFp = inv(Fp);
%L2a = (over_outer(iFp,iFp)*L2e*iJ11 + 2*(over_outer(iFp_old, iFp*S2b) + ...
%    over_outer(iFp*S2b, iFp_old)))*iJ21*over_outer(iFp',iFp');

end

% This function is used for the numerical jacobian
function [R,Phi,Phid] = Res2(mp,theta,Dt,rhop0,k_old,rhop_old,eb_tr_e,eb_tr_vol,X)
    tau_y  = mp.tau_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
    t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
    H      = mp.Hd     *exptemp_fun(mp.gammaH,theta,mp.thetaH);
    f      = 1+mp.alpha*(theta-mp.theta0);

    rhop = X(1); Dy = X(2); eb_e = X(3); eb_vol = X(4);
    ebb_vol = eb_vol - 3*log(f);

    q = rhop0/((1-rhop)*rhop^2)^(1/3);
    Meb = eb_e*(3*mp.G);
    Mmb = ebb_vol*mp.Kb;
    Mmbh = Mmb - mp.B*q;

    a1 = 1+(mp.c1-1)*(1-rhop)^mp.n1;
    a2 = mp.c2*(1-rhop)^mp.n2;

    Phi  = 1/(3*mp.G) * (a1*Meb^2  + a2*Mmbh^2 - tau_y^2);
    %Phid = 1/(3*mp.G) * (a1*Meb^2 + a2*Mmbh^2 - mp.tau_yd^2);
    Phid = -1/(3*mp.G)*mp.tau_yd^2;
    if Phid > 0
        warning('Inaccessible state!')
    %    keyboard
    end
    eta = 0;
    if Phi > 0 
        eta = mp.eta0*(-Phi/Phid)^mp.n;
    end

    % Residuals
    % Derivatives of the yield surfaces.
    dPdMeb = 2*a1*eb_e;
    dPdMmbh = 2*a2*Mmbh/(3*mp.G);

    % Residuals
    Rp = log(rhop) - log(rhop_old) + Dy*dPdMmbh;
    Rg = t_star*Dy - Dt*eta;
    Re = eb_e - eb_tr_e + Dy*dPdMeb;
    Rm = eb_vol - eb_tr_vol + Dy*dPdMmbh;
    R = [Rp;Rg;Re;Rm];
end

% For temp. dependent material properties
function e = exptemp_fun(gamma,theta,thetaX)
    if theta < thetaX
        e = 1;
    else
        e = exp(gamma/theta*(1-theta/thetaX));
    end
end

