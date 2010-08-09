function d=dblock0(p,x1,x2,y1,y2,z1,z2)
dx1 =  x1-p(:,1);
dx2 = -x2+p(:,1);
dy1 =  y1-p(:,2);
dy2 = -y2+p(:,2);
dz1 =  z1-p(:,3);
dz2 = -z2+p(:,3);

% Volume
d = -min(min(min(min(min(-dx1,-dx2),-dy1),-dy2),-dz1),-dz2);
%d = max([dx1,dx2,dy1,dy2,dz1,dz2],2);

% Edges
%%{
e1  = sqrt(dx1.^2+dy1.^2);
e2  = sqrt(dx1.^2+dy2.^2);
e3  = sqrt(dx2.^2+dy1.^2);
e4  = sqrt(dx2.^2+dy2.^2);
e5  = sqrt(dx1.^2+dz1.^2);
e6  = sqrt(dx1.^2+dz2.^2);
e7  = sqrt(dx2.^2+dz1.^2);
e8  = sqrt(dx2.^2+dz2.^2);
e9  = sqrt(dy1.^2+dz1.^2);
e10 = sqrt(dy1.^2+dz2.^2);
e11 = sqrt(dy2.^2+dz1.^2);
e12 = sqrt(dy2.^2+dz2.^2);
ix = dx1>0 & dy1>0; d(ix) = e1(ix);
ix = dx1>0 & dy2>0; d(ix) = e2(ix);
ix = dx2>0 & dy1>0; d(ix) = e3(ix);
ix = dx2>0 & dy2>0; d(ix) = e4(ix);
ix = dx1>0 & dz1>0; d(ix) = e5(ix);
ix = dx1>0 & dz2>0; d(ix) = e6(ix);
ix = dx2>0 & dz1>0; d(ix) = e7(ix);
ix = dx2>0 & dz2>0; d(ix) = e8(ix);
ix = dy1>0 & dz1>0; d(ix) = e9(ix);
ix = dy1>0 & dz2>0; d(ix) = e10(ix);
ix = dy2>0 & dz1>0; d(ix) = e11(ix);
ix = dy2>0 & dz2>0; d(ix) = e12(ix);
%}

% Corners
%%{
d1 = sqrt(dx1.^2+dy1.^2+dz1.^2);
d2 = sqrt(dx2.^2+dy1.^2+dz1.^2);
d3 = sqrt(dx1.^2+dy2.^2+dz1.^2);
d4 = sqrt(dx2.^2+dy2.^2+dz1.^2);
d5 = sqrt(dx1.^2+dy1.^2+dz2.^2);
d6 = sqrt(dx2.^2+dy1.^2+dz2.^2);
d7 = sqrt(dx1.^2+dy2.^2+dz2.^2);
d8 = sqrt(dx2.^2+dy2.^2+dz2.^2);
ix = dx1>0 & dy1>0 & dz1>0; d(ix) = d1(ix);
ix = dx2>0 & dy1>0 & dz1>0; d(ix) = d2(ix);
ix = dx1>0 & dy2>0 & dz1>0; d(ix) = d3(ix);
ix = dx2>0 & dy2>0 & dz1>0; d(ix) = d4(ix);
ix = dx1>0 & dy1>0 & dz2>0; d(ix) = d5(ix);
ix = dx2>0 & dy1>0 & dz2>0; d(ix) = d6(ix);
ix = dx1>0 & dy2>0 & dz2>0; d(ix) = d7(ix);
ix = dx2>0 & dy2>0 & dz2>0; d(ix) = d8(ix);
%}


