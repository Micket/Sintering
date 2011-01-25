function d=drectangle0(p,x1,x2,y1,y2)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

dx1 =  x1-p(:,1);
dx2 = -x2+p(:,1);
dy1 =  y1-p(:,2);
dy2 = -y2+p(:,2);

d1 = sqrt(dy1.^2+dx1.^2);
d2 = sqrt(dy1.^2+dx2.^2);
d3 = sqrt(dy2.^2+dx1.^2);
d4 = sqrt(dy2.^2+dx2.^2);

d=max([dx1,dx2,dy1,dy2],[],2);
ix = dy1>0 & dx1>0; d(ix) = d1(ix);
ix = dy1>0 & dx2>0; d(ix) = d2(ix);
ix = dy2>0 & dx1>0; d(ix) = d3(ix);
ix = dy2>0 & dx2>0; d(ix) = d4(ix);
