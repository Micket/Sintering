function [S,state,La,conv,debug] = sintering_matmod(C_old,C,Dt,theta,state,mp,print,numerical)
%function [S,state,La,conv] = sintering_matmod(C_old,C,Dt,theta,state,mp,print)
d = size(C,1);
I = eye(d);
Idev = over_outer(I,I) - outer(1/d*I,I);

maxiter = 20;
tolerance = [1e-12,1e-10,1e-12,1e-12]';

tau_y  = mp.tau_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
f      = 1+mp.alpha*(theta-mp.theta0);

state_old = state;
rhop0     = state.rhop0;

Cb_old = (state_old.Fp')\C_old/state_old.Fp;
Eb_old = 1/2*logm(Cb_old); %1/2*(Cb_old - I);
Eb_old_vol = trace(Eb_old);
Eb_old_dev = Eb_old - Eb_old_vol/d*I;
Eb_old_e = sqrt(2/3)*norm(Eb_old_dev);

Cb_tr = (state_old.Fp')\C/state_old.Fp;
Eb_tr = 1/2*logm(Cb_tr); %1/2*(Cb_tr - I);
Eb_tr_vol = trace(Eb_tr);
Eb_tr_dev = Eb_tr - Eb_tr_vol/d*I;
Eb_tr_e = sqrt(2/3)*norm(Eb_tr_dev);

X = [state_old.rhop,0,Eb_old_e,Eb_old_vol]';

for j = 1:maxiter
    % Analytical jacobian
    rhop = X(1); Dy = X(2); Eb_e = X(3); Eb_vol = X(4);
    Ebb_vol = Eb_vol - d*log(f);
    
    % Common expressions (and derivatives w.r.t. rhop)
    q = rhop0/((1-rhop)*rhop^2)^(1/3); qp  = rhop0*(3*rhop-2)/(3*(1-rhop)^(4/3)*rhop^(5/3));
    a1 = 1+(mp.c1-1)*(1-rhop)^mp.n1;   a1p = -(mp.c1-1)*mp.n1*(1-rhop)^(mp.n1-1);
    a2 = mp.c2*(1-rhop)^mp.n2;         a2p = -mp.c2*mp.n2*(1-rhop)^(mp.n2-1);

    Meb = Eb_e*3*mp.G;
    Mmb = Ebb_vol*mp.Kb;
    Mmbh = Mmb - mp.B*q;

    Phi  = 1/mp.const * (a1*Meb^2  + a2*Mmbh^2 - tau_y^2);
    eta = 0;
    etap = 0;
    if Phi > 0 
        eta = (Phi/mp.tau_c)^mp.n;
        etap = mp.n*(Phi/mp.tau_c)^(mp.n-1);
    end

    % Derivatives of the yield surface
    dPdrhop = 1/mp.const * (a1p*Meb^2   + (a2p*Mmbh^2 - 2*a2*Mmbh*mp.B*qp));
    dPdMeb  = 2*a1*Meb/mp.const;
    dPdMmbh = 2*a2*Mmbh/mp.const;
    dPdMmb  = dPdMmbh;

    % Residuals
    Rp = log(rhop) - log(state_old.rhop) + Dy*dPdMmbh;
    Rg = t_star*Dy - Dt*eta;
    Re = Eb_e - Eb_tr_e + Dy*dPdMeb;
    Rm = Eb_vol - Eb_tr_vol + Dy*dPdMmbh;
    R = [Rp;Rg;Re;Rm];

    if print
        fprintf('Constitutive iteration %2d, Res = ',j); fprintf('%11.4e, ',R); fprintf('\n');
        %fprintf('                           X   = ',j); fprintf('%11.4e, ',X); fprintf('\n');
    end
    % Jacobian
    Jp = [1/rhop + 2*Dy/mp.const*(a2p*Mmbh-a2*mp.B*qp),...
          dPdMmbh,...
          0,...
          Dy*2*a2*mp.Kb/mp.const];

    h = -Dt*etap/mp.tau_c;
    Jg = [h*dPdrhop,...
          t_star,...
          h*(3*mp.G)*dPdMeb,...
          h*mp.Kb*dPdMmb];

    Je = [2*Dy*a1p*Eb_e,...
          dPdMeb,...
          1+2*Dy*a1*3*mp.G/mp.const,...
          0];

    Jm = [2*Dy/mp.const*(a2p*Mmbh-a2*mp.B*qp),...
          dPdMmbh,...
          0,...
          1+Dy*2*a2*mp.Kb/(3*mp.G)];
          
    J = [Jp; Jg; Je; Jm];

    conv = all(abs(R) <= tolerance);
    if conv
        break
    end

    % Debugging; numerical jacobian
    %%{
    Rn = Res2(mp,theta,Dt,state_old,Eb_tr_e,Eb_tr_vol,X);
    hn = 1e-6*[1e-6,1e-6,1e-6,1e-6]';
    Jn = zeros(length(hn));
    for i = 1:length(hn);
        Xh = X; Xh(i) = Xh(i) + hn(i);
        Rh = Res2(mp,theta,Dt,state_old,Eb_tr_e,Eb_tr_vol,Xh);
        Jn(:,i) = (Rh-Rn)/hn(i);
    end
    % Check for any difference
    dJ = Jn-J;
    dR = Rn-R;
    if norm(dJ(:))/norm(J(:)) > 1e-0 || norm(dR) > 1e-1
        disp('Numerical jacobian doesn''t match analytical!');
        dJ
        keyboard
    end
    %}

    X = X -J\R;
end

if print
    fprintf('Constitutive iterations = %2d, Dy = %e\n',j,Dy)
end
if ~conv
    warning(['Material not converged, res = ',num2str(R')])
end
% Postprocessing, calculating the full forms
if Eb_e == 0
    Eb_dev = zeros(d);
else
    %Eb_dev = Eb_e/(Eb_e + Dy*dPdMeb)*Eb_tr_dev;
    Eb_dev = Eb_e/Eb_tr_e*Eb_tr_dev;
end
Mb = 2*mp.G*Eb_dev + Mmb*I;
nu = 2*a1*Eb_dev + 1/d*dPdMmbh*I;
iAb = expm(-Dy*nu);
Fp = iAb\state_old.Fp;
Cb = (Fp')\C/Fp;
Cb2 = expm(2*(Eb_dev + Eb_vol/d*I));
%Cb2 = 2*(Eb_dev + Eb_vol/d*I) + I;
if norm(Cb - Cb2) > 1e-5
    disp('Inconsistencies in Cb!')
    Cb - Cb2
    keyboard
end

Sb = Mb; % it's actually Cb\Mb, but Mb has been approximated as Sb
S = Fp\Sb/(Fp');

state.Fp   = Fp;
state.rhop = rhop;

%%%%%%%%%%%%%%%%%%% ATS Tensor %%%%%%%%%%%%%%%%%%%%%%
I4 = eye(d*d);
Iv = zeros(d*d,1); Iv(1:d) = 1;
% Elastic part
Le = 2*mp.G*Idev + mp.Kb*Iv*Iv';
% Plastic part
iFp = inv(Fp);
iFp_old = inv(state_old.Fp);
dCbtrdC = over_outer(iFp_old',iFp_old');
dEbtrvoldC = 1/2*dCbtrdC*m2v(I);
if Eb_tr_e == 0
    dEbtredC = 1/3*dCbtrdC*m2v(1-I)*sqrt(2/3); % TODO CHECK THIS! VERY UNSURE ON HOW TO HANDLE IT
    %dEbtredC = 1/3*dCbtrdC*m2v(zeros(size(I))); % I *think* it should not matter
else
    dEbtredC = 1/(3*Eb_tr_e)*dCbtrdC*m2v(Eb_tr_dev);
end
dRdC = [zeros(1,d*d);zeros(1,d*d);-dEbtredC';-dEbtrvoldC'];
D = -J\dRdC;
D_rho = v2m(D(1,:));
D_gamma = v2m(D(2,:));
D_e = v2m(D(3,:));
D_m = v2m(D(4,:));

dEbtrdevdC = 1/2*Idev*dCbtrdC;
dEbdevdC = 2*(1+2*a1*Dy)^(-2)*outer(Eb_tr_dev,a1*D_gamma+a1p*Dy*D_rho) + 1/(1+2*a1*Dy)*dEbtrdevdC;
2*a1*(1+2*a1*Dy)^(-2)*outer(Eb_tr_dev,D_gamma) + 1/(1+2*a1*Dy)*dEbtrdevdC;
dEbvoldC = outer(1/d*I,D_m);
dCbdC = 2*(dEbdevdC + dEbvoldC);

dMmbhdC = mp.Kb*D_m - mp.B*qp*D_rho;
dnudC = 2*(outer(a1p*Eb_dev,D_rho) + a1*dEbdevdC + 1/(3*mp.G*d)*outer(I,(a2p*Mmbh*D_rho+a2*dMmbhdC)));
diAbdC = -over_outer(iAb,I)*(outer(nu,D_gamma) + Dy*dnudC);

La = 1/2*over_outer(iFp,iFp)*Le*dCbdC + ...
    (over_outer(iFp_old, iFp*Sb) + over_outer(iFp*Sb, iFp_old))*diAbdC;

% Numerical ATS
debug.Cb     = Cb;
debug.Cb_tr  = Cb_tr;
debug.Eb_tr  = Eb_tr;
debug.Eb_dev = Eb_dev;
debug.Eb_vol = Eb_vol;
debug.Eb_tr_dev = Eb_tr_dev;
debug.Eb_tr_vol = Eb_tr_vol;
debug.iAb    = iAb;
debug.nu     = nu;
debug.rhop   = rhop;
debug.Eb_vol = Eb_vol;
debug.Eb_e   = Eb_e;
debug.Dy     = Dy;
debug.iterations = j;

%%{
if nargin <= 7
    h = 1e-11;
    offdiag = (d*d-d)/2;
    Lan      = zeros(d*d);
    dCbdCn   = zeros(d*d);
    dCbtrdCn = zeros(d*d);
    dEbdevdCn = zeros(d*d);
    dEbvoldCn = zeros(d*d);
    dEbtrdCn = zeros(d*d);
    dEbtrdevdCn = zeros(d*d);
    diAbdCn  = zeros(d*d);
    dnudCn   = zeros(d*d);
    D_mn     = zeros(d*d,1);
    D_en     = zeros(d*d,1);
    D_rhon   = zeros(d*d,1);
    D_gamman = zeros(d*d,1);
    for q = 1:d*d
        Ch = m2v(C);
        Ch(q) = Ch(q) + h;
        Ch = v2m(Ch);
        Ch = 1/2*(Ch+Ch');
        [Sh,stateh,Lah,convh,debugh] = sintering_matmod(C_old,Ch,Dt,theta,state_old,mp,false,false);
        if ~convh
            disp('numerical derivative didn''t converge');
            break
            %keyboard
        end
        Lan(:,q)       = m2v(Sh-S)/h;
        dCbdCn(:,q)    = m2v(debugh.Cb-Cb)/h;
        dCbtrdCn(:,q)  = m2v(debugh.Cb_tr-Cb_tr)/h;
        dEbdevdCn(:,q) = m2v(debugh.Eb_dev-Eb_dev)/h;
        dEbvoldCn(:,q) = m2v(I/d*(debugh.Eb_vol-Eb_vol))/h;
        dEbtrdCn(:,q)  = m2v(debugh.Eb_tr-Eb_tr)/h;
        dEbtrvoldCn(:,q) = m2v(I/d*(debugh.Eb_tr_vol-Eb_tr_vol))/h;
        dEbtrdevdCn(:,q) = m2v(debugh.Eb_tr_dev-Eb_tr_dev)/h;
        diAbdCn(:,q) = m2v(debugh.iAb-iAb)/h;
        dnudCn(:,q)  = m2v(debugh.nu-nu)/h;
        D_mn(q)      = m2v(debugh.Eb_vol - Eb_vol)/h;
        D_en(q)      = m2v(debugh.Eb_e - Eb_e)/h;
        D_rhon(q)    = m2v(debugh.rhop - rhop)/h;
        D_gamman(q)  = m2v(debugh.Dy - Dy)/h;
    end
    D_mn = v2m(D_mn);
    D_en = v2m(D_en);
    D_rhon = v2m(D_rhon);
    D_gamman = v2m(D_gamman);
    Lan(:,d+1:d+offdiag) = Lan(:,d+1:d+offdiag) + Lan(:,end+1-offdiag:end);
     
    %Lan(d+1:end,:) = Lan(d+1:end,:)*2;
    %Lan(:,d+1:end) = Lan(:,d+1:end)*2;
    
    err = (Lan(1:3,1:3) - La(1:3,1:3))/norm(Lan(1:3,1:3));
    if q==d*d && norm(err) > 0.1
        Lan
        La
        keyboard
    end
    %La = Lan;
end
%}

end

% This function is used for the numerical jacobian
%%{
function [R,Phi] = Res2(mp,theta,Dt,state_old,Eb_tr_e,Eb_tr_vol,X)
    tau_y  = mp.tau_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
    t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
    f      = 1+mp.alpha*(theta-mp.theta0);
    d = size(state_old.Fp,1);

    rhop = X(1); Dy = X(2); Eb_e = X(3); Eb_vol = X(4);
    Ebb_vol = Eb_vol - d*log(f);

    q = state_old.rhop0/((1-rhop)*rhop^2)^(1/3);
    Meb = Eb_e*(3*mp.G);
    Mmb = Ebb_vol*mp.Kb;
    Mmbh = Mmb - mp.B*q;

    a1 = 1+(mp.c1-1)*(1-rhop)^mp.n1;
    a2 = mp.c2*(1-rhop)^mp.n2;

    Phi  = 1/mp.const * (a1*Meb^2  + a2*Mmbh^2 - tau_y^2);
    eta = 0;
    if Phi > 0 
        eta = (Phi/mp.tau_c)^mp.n;
    end

    % Derivatives of the yield surfaces.
    dPdMeb  = 2*a1*Meb/mp.const;
    dPdMmbh = 2*a2*Mmbh/mp.const;

    % Residuals
    Rp = log(rhop) - log(state_old.rhop) + Dy*dPdMmbh;
    Rg = t_star*Dy - Dt*eta;
    Re = Eb_e - Eb_tr_e + Dy*dPdMeb;
    Rm = Eb_vol - Eb_tr_vol + Dy*dPdMmbh;
    R = [Rp;Rg;Re;Rm];
end
%}

% For temp. dependent material properties
function e = exptemp_fun(gamma,theta,thetaX)
    if theta < thetaX
        e = 1;
    else
        e = exp(gamma/theta*(1-theta/thetaX));
    end
end

