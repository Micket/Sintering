function h=huniform(p,varargin)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
h = interp1([-0.1,0.5,1.1]',[10,1,10]',p(:,1));
%h=ones(size(p,1),1);
