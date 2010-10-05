function [stress,elasticity] = yeoh(Cb,modpar,lagrangian)
mu = modpar.mu;
lambda = modpar.lambda;
c2 = modpar.c2;
c3 = modpar.c3;
d = size(Cb,1);
I = eye(d);
J = sqrt(det(Cb));

if lagrangian
    C = Cb;
    I_C = trace(C);
    iC = inv(C);
    J  = sqrt(det(C));
    I2 = (over_outer(iC,iC) + under_outer(iC,iC));
    stress = mu*(eye(d)-iC)+lambda*log(J)*iC + (4*c2*(I_C-d)+6*c3*(I_C-d)^2)*eye(d);
    elasticity = lambda*outer(iC,iC) + (mu-lambda*log(J))*I2 + (8*c2+24*c3*(I_C-d))*outer(I,I);
else
    b = Cb;
    I_b = trace(b);
    i2 = (over_outer(I,I) + under_outer(I,I));
    stress = mu/J*(b-I) + lambda/J*log(J)*I + 1/J*(4*c2*(I_b-d)+6*c3*(I_b-d)^2)*b;
    elasticity = lambda/J*outer(I,I) + (mu-lambda*log(J))/J*i2 + (8*c2+24*c3*(I_b-d))/J*outer(b,b);
end
end

