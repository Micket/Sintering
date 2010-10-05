function [S,state,Ea,conv,debug] = sintering_small(E_old,E,Dt,theta_old,theta,state,mp,print,numerical)
% Notation: S = stress (sigma), E = strain (epsilon)

d = size(E,1);
I = eye(d);
Idev = over_outer(I,I) - outer(1/d*I,I);
Iv = m2v(I);
IoI = Iv*Iv';

Ee = 2*mp.G*Idev + mp.Kb*IoI;

maxiter = 20;
tolerance = [1e-11,1e-9,1e-11,1e-11]';

S_y    = mp.S_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
f      = mp.alpha*(theta-mp.theta0);
f_old  = mp.alpha*(theta_old-mp.theta0);

state_old = state;

Deps_free = (f-f_old)*I;
Eb_0 = E_old - state_old.E_p + Deps_free;
Eb_tr = E - state_old.E_p;

Eb_0_vol = trace(Eb_0);
Eb_0_dev = Eb_0 - Eb_0_vol/d*I;
Eb_0_e = sqrt(2/3)*norm(Eb_0_dev);

Eb_tr_vol = trace(Eb_tr);
Eb_tr_dev = Eb_tr - Eb_tr_vol/d*I;
Eb_tr_e = sqrt(2/3)*norm(Eb_tr_dev);

X = [state_old.rhop,0,Eb_0_e,Eb_0_vol]';
% Not a good guess, the trail values;
%X = [state_old.rhop,0,Eb_tr_e,Eb_tr_vol]';

for j = 1:maxiter
    % Analytical jacobian
    rhop = X(1); Dy = X(2); Eb_e = X(3); Eb_vol = X(4);
    
    % Common expressions (and derivatives w.r.t. rhop)
    %q = state.rhop0/((1-rhop)*rhop^2)^(1/3);
    %qp = state.rhop0*(3*rhop-2)/(3*(1-rhop)^(4/3)*rhop^(5/3));
    q = (rhop/(1-rhop))^(1/3);
    qp = 1/(3*(rhop/(1-rhop))^(2/3)*(rhop-1)^2);
    a = mp.c*(1-rhop);
    ap = -mp.c;

    Se = Eb_e*3*mp.G;
    Sm = (Eb_vol - d*f)*mp.Kb;
    Smh = Sm - mp.B*q;

    Phi  = 1/mp.P * (Se^2  + a*Smh^2 - S_y^2);
    eta = 0;
    etap = 0;
    if Phi > 0 
        eta = (Phi/mp.S_c)^mp.n;
        etap = mp.n/mp.S_c*(Phi/mp.S_c)^(mp.n-1);
    end

    % Derivatives of the yield surface
    dPdrhop = 1/mp.P * (Se^2   + (ap*Smh^2 - 2*a*Smh*mp.B*qp));
    dPdSe  = 2*Se/mp.P;
    dPdSmh = 2*a*Smh/mp.P;
    dPdSm  = dPdSmh;

    % Residuals
    Rp = log(rhop) - log(state_old.rhop) + Dy*dPdSmh;
    Rg = t_star*Dy - Dt*eta;
    Re = Eb_e - Eb_tr_e + Dy*dPdSe;
    Rm = Eb_vol - Eb_tr_vol + Dy*dPdSmh;
    R = [Rp;Rg;Re;Rm];

    if print
        fprintf('Constitutive iteration %2d, Res = ',j); fprintf('%11.4e, ',R); fprintf('\n');
        %fprintf('                           X   = ',j); fprintf('%11.4e, ',X); fprintf('\n');
    end
    % Jacobian
    Jp = [1/rhop + 2*Dy/mp.P*(ap*Smh - a*mp.B*qp),...
          dPdSmh,...
          0,...
          Dy*2*a*mp.Kb/mp.P];

    Jg = [-Dt*etap*dPdrhop,...
          t_star,...
          -Dt*etap*(3*mp.G)*dPdSe,...
          -Dt*etap*mp.Kb*dPdSm];

    Je = [0,...
          dPdSe,...
          1+2*Dy*3*mp.G/mp.P,...
          0];

    Jm = [2*Dy/mp.P*(ap*Smh - a*mp.B*qp),...
          dPdSmh,...
          0,...
          1+Dy*2*a*mp.Kb/mp.P];
          
    J = [Jp; Jg; Je; Jm];

    conv = all(abs(R) <= tolerance);
    if conv
        break
    end

    % Debugging; numerical jacobian
    %{
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
    J = Jn;
    %}

    X = X -J\R;
end
conv = all(abs(R) <= tolerance);
if ~conv
    warning(['Material not converged, res = ',num2str(abs(R)' <= tolerance')])
    keyboard
end
% Postprocessing, calculating the full forms
F = (1+2*Dy*3*mp.G/mp.P);
if Eb_e == 0
    Eb_dev = zeros(d);
else
    Eb_dev = 1/F*Eb_tr_dev;
end
S = 2*mp.G*Eb_dev + Sm*I;
nu = 2*Eb_dev + 1/d*dPdSmh*I;
E_p = state_old.E_p + Dy*nu;

state.E_p  = E_p;
state.rhop = rhop;

%%%%%%%%%%%%%%%%%%% ATS Tensor %%%%%%%%%%%%%%%%%%%%%%
dEbtrvoldE = Iv;
if Eb_tr_e == 0
    dEbtredE = zeros(d*d,1); % I *think* it doesn't matter
else
    dEbtredE = 2/(3*Eb_tr_e)*m2v(Eb_tr_dev);
end
dRdE = [zeros(1,d*d);zeros(1,d*d);-dEbtredE';-dEbtrvoldE'];
D = -J\dRdE;
D_rho = v2m(D(1,:));
D_gamma = v2m(D(2,:));
D_e = v2m(D(3,:));
D_v = v2m(D(4,:));

if D_v(1,1) < 0.1
%    keyboard
end
%D_v
dEbvoldE = outer(I/d,D_v);
dEbdevdE = outer(2*3*mp.G/mp.P/F^2*Eb_tr_dev,D_gamma) + 1/F*Idev;
dEbdE = dEbdevdE + dEbvoldE;
Ea = Ee*dEbdE;

% Numerical ATS
debug.rhop   = rhop;
debug.Eb_vol = Eb_vol;
debug.Eb_e   = Eb_e;
debug.Dy     = Dy;

%{
if nargin <= 7
    offdiag = (d*d-d)/2;
    Ean      = zeros(d*d);
    D_mn     = zeros(d*d,1);
    D_en     = zeros(d*d,1);
    D_rhon   = zeros(d*d,1);
    D_gamman = zeros(d*d,1);
    h = 1e-10;
    for q = 1:d*d
        Eh    = m2v(E);
        Eh(q) = Eh(q) + h;
        Eh    = v2m(Eh);
        %Eh    = 1/2*(Eh+Eh');
        [Sh,stateh,Lah,convh,debugh] = sintering_small(E_old,Eh,Dt,theta,state_old,mp,false,false);
        if ~convh
            disp('Numerical ATS didn''t converge.');
            break
            %keyboard
        end
        Ean(:,q)     = m2v(Sh-S)/h;
        D_mn(q)      = m2v(debugh.Eb_vol - Eb_vol)/h;
        D_en(q)      = m2v(debugh.Eb_e - Eb_e)/h;
        D_rhon(q)    = m2v(debugh.rhop - rhop)/h;
        D_gamman(q)  = m2v(debugh.Dy - Dy)/h;
    end
    D_mn = v2m(D_mn);
    D_en = v2m(D_en);
    D_rhon = v2m(D_rhon);
    D_gamman = v2m(D_gamman);
    %Ean(:,d+1:d+offdiag) = Ean(:,d+1:d+offdiag) + Ean(:,end+1-offdiag:end);
     
    err = (Ean(1:3,1:3) - Ea(1:3,1:3))/norm(Ean(1:3,1:3));
    if q==d*d && norm(err) > 1
        Ean
        Ea
        keyboard
    end
    Ea = Ean;
end
%}

end

% This function is used for the numerical jacobian
%{
function [R,Phi] = Res2(mp,theta,Dt,state_old,Eb_tr_e,Eb_tr_vol,X)
    S_y    = mp.S_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
    t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
    f      = mp.alpha*(theta-mp.theta0);
    d = size(state_old.E_p,1);

    rhop = X(1); Dy = X(2); Eb_e = X(3); Eb_vol = X(4);

    q = state_old.rhop0/((1-rhop)*rhop^2)^(1/3);
    Se = Eb_e*(3*mp.G);
    Sm = (Eb_vol-d*f)*mp.Kb;
    Smh = Sm - mp.B*q;

    a = mp.c*(1-rhop);

    Phi  = 1/mp.P * (Se^2  + a*Smh^2 - S_y^2);
    eta = 0;
    if Phi > 0 
        eta = (Phi/mp.S_c)^mp.n;
    end

    % Derivatives of the yield surfaces.
    dPdSe  = 2*Se/mp.P;
    dPdSmh = 2*a*Smh/mp.P;

    % Residuals
    Rp = log(rhop) - log(state_old.rhop) + Dy*dPdSmh;
    Rg = t_star*Dy - Dt*eta;
    Re = Eb_e - Eb_tr_e + Dy*dPdSe;
    Rm = Eb_vol - Eb_tr_vol + Dy*dPdSmh;
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

