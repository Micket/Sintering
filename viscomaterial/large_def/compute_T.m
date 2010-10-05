function [Ta] = compute_T(x,t,dN,sigma)
dv = t*abs(x(1,2)*x(2,1)-x(1,1)*x(2,2)-x(1,2)*x(3,1)+...
    x(2,2)*x(3,1)+x(1,1)*x(3,2)-x(2,1)*x(3,2))/2;
Ta = dv*reshape(sigma*dN,[],1);
