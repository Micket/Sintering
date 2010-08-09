function h = h_4particle_rve(p,x,r)

a0 = sqrt(r^2-x^2);

dc1 = dcircle(p,-x,-x,r);
dc2 = dcircle(p, x,-x,r);
dc3 = dcircle(p, x, x,r);
dc4 = dcircle(p,-x, x,r);

d = min([dc1,dc2,dc3,dc4],[],2); % Union of dc.

% The points, at a0
y = [0,-x+a0;...
     -x+a0,0;...
     0,x-a0;...
     x-a0,0];

% Right
dp = -sqrt( (p(:,1)-x+a0).^2 + p(:,2).^2 );
dm = min(d(:,1),dp); 
ii = p(:,2)*a0 < (p(:,1)-x+a0)*x & p(:,2)*a0 > -(p(:,1)-x+a0)*x; d(ii) = dm(ii);

% Left
dp = -sqrt( (p(:,1)+x-a0).^2 + p(:,2).^2 );
dm = min(d(:,1),dp); 
ii = p(:,2)*a0 < -(p(:,1)+x-a0)*x & p(:,2)*a0 > (p(:,1)+x-a0)*x; d(ii) = dm(ii);

% Bottom
dp = -sqrt( p(:,1).^2 + (p(:,2)-x+a0).^2 );
dm = min(d(:,1),dp); 
ii = p(:,1)*a0 < (p(:,2)-x+a0)*x & p(:,1)*a0 > -(p(:,2)-x+a0)*x; d(ii) = dm(ii);

% Top
dp = -sqrt( p(:,1).^2 + (p(:,2)+x-a0).^2 );
dm = min(d(:,1),dp);
ii = p(:,1)*a0 < (-p(:,2)-x+a0)*x & p(:,1)*a0 > -(-p(:,2)-x+a0)*x; d(ii) = dm(ii);

%dp = sqrt( (p(:,1)-y(1,1)).^2 + (p(:,2)-y(1,2)).^2 );
%dm = max(dc(:,1),dp); 
%ii = p(:,1) > x-a0; d(ii) = dm(ii);

%h = 1./(1+d);

h = blä...
