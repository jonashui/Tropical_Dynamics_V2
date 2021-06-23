clear
close all
b = 0;
F = [0 -1 1; 0 1 0];
G = [b-1 0 -1; -1 1 -1];
dF = [1 -1];
dG = [1 -1];
ulim = 3*[-1 1];
vlim = -1+ulim;
w0 = [-2 -3];
TropicalPhasePlane(F,G,dF,dG,ulim,vlim)
EpsilonPhasePlot(F,G,dF,dG,w0)