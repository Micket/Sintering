function d=dcirc_special(p,xc,yc,r)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.


rp =sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2);

d = rp - r;
d(rp<0.75*r) = 0.5*r - rp(rp<0.75*r);
d(rp<0.50*r) = rp(rp<0.5*r) - 0.5*r;
