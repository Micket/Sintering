function [stress,elasticity] = neo_hook(F,modpar,lagrangian)
mu = modpar.mu;
lambda = modpar.lambda;
d = size(F,1);
I = eye(d);
J = det(F);
if lagrangian
    C = F'*F;
    iC = inv(C);
    stress = mu*(eye(d)-iC)+lambda*log(J)*iC;
    I2 = over_outer(iC,iC) + under_outer(iC,iC);
    elasticity = lambda*outer(iC,iC) + (mu-lambda*log(J))*I2;
else
    b = F*F';
    stress = mu/J*(b-I)+lambda/J*log(J)*I;
    i2 = over_outer(I,I)+under_outer(I,I);
    elasticity = lambda/J*outer(I,I) + (mu-lambda*log(J))/J*i2;
end
end

