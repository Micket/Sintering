function [Kc,Ks] = compute_K(dN,x,t,c,sigma)
dv = t*abs(x(1,2)*x(2,1)-x(1,1)*x(2,2)-x(1,2)*x(3,1)+x(2,2)*x(3,1)+x(1,1)*x(3,2)-x(2,1)*x(3,2))/2;
[d,s] = size(dN);
Kc = zeros(d*s); Ks = zeros(d*s);
if d == 1
    q = [1];
elseif d == 2
    q = [1,3;4,2];
else
    q = [1,4,6;7,2,5;9,8,3];
end
for i = 1:d; for j = 1:d
    Kc(i:d:end,j:d:end) = dv*dN'*c(q(i,:),q(j,:))*dN;
end; end
Ks_tmp = dN'*sigma*dN*dv;
for ij = 1:d;
    Ks(ij:d:end,ij:d:end) = Ks_tmp;
end

