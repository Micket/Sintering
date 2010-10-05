function [S,state,La,conv,debug] = sintering_matmod(C_old,C,Dt,theta,state,mp,print,numerical)
%function [S,state,La,conv] = sintering_matmod(C_old,C,Dt,theta,state,mp,print)
d = size(C,1);
I = eye(d);
I4 = eye(d*d);
Iv = m2v(I);
IoI = Iv*Iv';
Idev = I4 - IoI/d;
% Elastic part
Le = 2*mp.G*Idev + mp.Kb*IoI;

maxiter = 20;
tolerance = 0;%1e-12;

tau_y  = mp.tau_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
f      = 1+mp.alpha*(theta-mp.theta0);

state_old = state;

Cb_old = state_old.iFp'*C_old*state_old.iFp;
Cb_tr = state_old.iFp'*C*state_old.iFp;

X = [m2v(Cb_old);Iv;state_old.rhop;0];

for j = 1:maxiter
    % Analytical jacobian
    Cb = v2m(X(1:d*d));
    A = v2m(X(1+d*d:d*d*2));
    rhop = X(end-1);
    Dy = X(end);

    Ebb = 1/2*(Cb/f^2-I);

    Ebb_vol = trace(Ebb);
    Ebb_dev = Ebb - Ebb_vol/d*I;
    Ebb_e = sqrt(2/3)*norm(Ebb_dev);
    
    % Common expressions (and derivatives w.r.t. rhop)
    q = state.rhop0/((1-rhop)*rhop^2)^(1/3); 
    qp = state.rhop0*(3*rhop-2)/(3*(1-rhop)^(4/3)*rhop^(5/3));
    a = mp.c*(1-rhop); 
    ap = -mp.c;

    % Stresses
    Meb = Ebb_e*3*mp.G;
    Mmb = Ebb_vol*mp.Kb;
    Mmbh = Mmb - mp.B*q;

    % Yield surface and overstress function
    Phi  = 1/mp.P * (Meb^2  + a*Mmbh^2 - tau_y^2);
    eta = 0;
    etap = 0;
    if Phi > 0 
        eta = (Phi/mp.tau_c)^mp.n;
        etap = mp.n/mp.tau_c*(Phi/mp.tau_c)^(mp.n-1);
    end

    % Derivatives of the yield surface
    dPdrhop = 1/mp.P * (Meb^2   + (ap*Mmbh^2 - 2*a*Mmbh*mp.B*qp));
    dPdMeb  = 2*Meb/mp.P;
    dPdMmbh = 2*a*Mmbh/mp.P;

    % Direction
    nu = 3/mp.P*2*mp.G*Ebb_dev + 1/d*dPdMmbh*I;

    % Residuals
    AtCbtr = A'*Cb_tr;
    Rc = Cb - AtCbtr*A;
    Ra = A - (I-Dy*nu);
    Rp = log(rhop) - log(state_old.rhop) + Dy*dPdMmbh;
    Rg = t_star*Dy - Dt*eta;
    R = [m2v(Rc);m2v(Ra);Rp;Rg];

    if print
        fprintf('Constitutive iteration %2d, Res = ',j); fprintf('%11.4e, ',R); fprintf('\n');
        %fprintf('                           X   = ',j); fprintf('%11.4e, ',X); fprintf('\n');
    end
    % Jacobian
    Jc = [I4,...
          -under_outer(I,AtCbtr) - over_outer(AtCbtr,I),...
          zeros(d*d,1),...
          zeros(d*d,1)];

    Ja = [Dy/mp.P*(3*mp.G*Idev + a*mp.Kb/d*IoI),...
          I4,...
          Dy/(mp.P*d)*(ap*Mmbh-a*mp.B*qp)*Iv,...
          m2v(nu)];

    Jp = [Dy*a*mp.Kb/mp.P*Iv',...
          zeros(1,d*d),...
          1/rhop + 2*Dy/mp.P*(ap*Mmbh - a*mp.B*qp),...
          dPdMmbh];

    Jg = [-Dt*etap/2*m2v(nu)'*Le,...
          zeros(1,d*d),...
          -Dt*etap*dPdrhop,...
          t_star];
          
    J = [Jc; Ja; Jp; Jg];

    conv = all(abs(R) <= tolerance);
    if conv
        break
    end

    % Debugging; numerical jacobian
    %%{
    hn = 1e-12;
    Jn = zeros(length(X));
    for i = 1:length(X);
        Xh = X; Xh(i) = Xh(i) + hn;
        Rh = Res2(mp,theta,Dt,state_old,Cb_tr,Xh);
        Jn(:,i) = (Rh-R)/hn;
    end
    % Check for any difference
    dJ = Jn-J;
    if norm(dJ(:))/norm(J(:)) > 1e-1
        disp('Numerical jacobian doesn''t match analytical!');
        dJ
        keyboard
    end
    %}

    X = X -J\R;
end
%if print
    fprintf('Constitutive iterations = %2d, Dy = %e\n',j,Dy)
%end
if ~conv
    warning(['Material not converged, res = ',num2str(R')])
end
iFp = state_old.iFp*A;
Cb = iFp'*C*iFp;

Mb = 2*mp.G*Ebb_dev + mp.Kb*Ebb_vol*I;
Sb = Mb; % Sb = Cb\Mb
S = iFp*Sb*(iFp');

state.iFp   = iFp;
state.rhop = rhop;

%%%%%%%%%%%%%%%%%%% ATS Tensor %%%%%%%%%%%%%%%%%%%%%%
% Plastic part
iFp_old = state_old.iFp;
dCbtrdC = over_outer(iFp_old',iFp_old');
dEbtrvoldC = 1/2*dCbtrdC*m2v(I);
iJ = inv(J);
Dc = iJ(1:d*d,1:d*d)*dCbtrdC;
Da = iJ(1+d*d:d*d*2,1:d*d)*dCbtrdC;
%Dr = ; % Can be calculated, but not needed.
%Dg = ;

La = 1/2*over_outer(iFp,iFp)*Le*Dc + ...
    (over_outer(iFp_old, iFp*Sb) + over_outer(iFp*Sb, iFp_old))*Da;

% Numerical ATS
debug.Cb = Cb;
debug.nu = nu;
debug.A = A;
debug.rhop = rhop;
debug.Dy = Dy;

%{
if nargin <= 7
    h = 1e-10;
    offdiag= (d*d-d)/2;
    Lan    = zeros(d*d);
    dnudCn = zeros(d*d);
    Dcn    = zeros(d*d);
    Dan    = zeros(d*d);
    Dpn    = zeros(1,d*d);
    Dgn    = zeros(1,d*d);
    for q = 1:d*d
        Ch = m2v(C);
        Ch(q) = Ch(q) + h;
        Ch = v2m(Ch);
        Ch = 1/2*(Ch+Ch');
        [Sh,stateh,Lah,convh,debugh] = sintering_full(C_old,Ch,Dt,theta,state_old,mp,false,false);
        if ~convh
            disp('numerical derivative didn''t converge');
            break
            %keyboard
        end
        Lan(:,q) = m2v(Sh-S)/h;
        dnudCn(:,q) = m2v(debugh.nu-nu)/h;
        Dcn(:,q) = m2v(debugh.Cb - Cb)/h;
        Dan(:,q) = m2v(debugh.A - A)/h;
        Dpn(q) = m2v(debugh.rhop - rhop)/h;
        Dgn(q) = m2v(debugh.Dy - Dy)/h;
    end
    Dpn = v2m(Dpn);
    Dgn = v2m(Dgn);
    Lan(:,d+1:d+offdiag) = Lan(:,d+1:d+offdiag) + Lan(:,end+1-offdiag:end);
     
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
function [R] = Res2(mp,theta,Dt,state_old,Cb_tr,X)
    tau_y  = mp.tau_yd *exptemp_fun(mp.gammaT,theta,mp.thetaT);
    t_star = mp.t_stard*exptemp_fun(mp.gammat,theta,mp.thetat);
    f      = 1+mp.alpha*(theta-mp.theta0);
    d = size(state_old.iFp,1);
    I = eye(d);

    Cb = v2m(X(1:d*d));
    A = v2m(X(1+d*d:d*d*2));
    rhop = X(end-1);
    Dy = X(end);

    Ebb = 1/2*(Cb/f^2-I);

    Ebb_vol = trace(Ebb);
    Ebb_dev = Ebb - Ebb_vol/d*I;
    Ebb_e = sqrt(2/3)*norm(Ebb_dev);
    
    % Common expressions (and derivatives w.r.t. rhop)
    q = state_old.rhop0/((1-rhop)*rhop^2)^(1/3); 
    qp = state_old.rhop0*(3*rhop-2)/(3*(1-rhop)^(4/3)*rhop^(5/3));
    a = mp.c*(1-rhop); 
    ap = -mp.c;

    % Stresses
    Meb = Ebb_e*3*mp.G;
    Mmb = Ebb_vol*mp.Kb;
    Mmbh = Mmb - mp.B*q;

    % Yield surface and overstress function
    Phi  = 1/mp.P * (Meb^2  + a*Mmbh^2 - tau_y^2);
    eta = 0;
    etap = 0;
    if Phi > 0 
        eta = (Phi/mp.tau_c)^mp.n;
        etap = mp.n/mp.tau_c*(Phi/mp.tau_c)^(mp.n-1);
    end

    % Derivatives of the yield surface
    dPdrhop = 1/mp.P * (Meb^2   + (ap*Mmbh^2 - 2*a*Mmbh*mp.B*qp));
    dPdMeb  = 2*Meb/mp.P;
    dPdMmbh = 2*a*Mmbh/mp.P;

    % Direction
    nu = 3/mp.P*2*mp.G*Ebb_dev + 1/d*dPdMmbh*I;

    % Residuals
    AtCbtr = A'*Cb_tr;
    Rc = Cb - AtCbtr*A;
    Ra = A - (I-Dy*nu);
    Rp = log(rhop) - log(state_old.rhop) + Dy*dPdMmbh;
    Rg = t_star*Dy - Dt*eta;
    R = [m2v(Rc);m2v(Ra);Rp;Rg];
end

% For temp. dependent material properties
function e = exptemp_fun(gamma,theta,thetaX)
    if theta < thetaX
        e = 1;
    else
        e = exp(gamma/theta*(1-theta/thetaX));
    end
end

