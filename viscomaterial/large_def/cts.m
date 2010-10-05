function [d0N,dN] = cts(X,x)
dsN = [-1,1,0;-1,0,1];
dsX = dsN*X;
dsx = dsN*x;
d0N = dsX\dsN;
dN  = dsx\dsN;
