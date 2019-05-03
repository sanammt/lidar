%%%%%%%% Start of Copyright Notice %%%%%%%
%
% copyright © 2018
% Written by Sanam Moslemi-Tabrizi
% All rights reserved
% Can NOT be copied or distrutated without including this copyright notice
%
%%%%%%%% End of Copyright Notice %%%%%%%%%%%
clear all
close all
syms x y z
a=1.55;
f=a/(x*sqrt(1-(y-a/z)^2))
gradient(f, [x, y, z])
gx=gradient(f, x)
gy=gradient(f, y)
gz=gradient(f, z)

Y1=2.6,
Z1=0.6
X1=10:5:100;
Neff_cap= a/Z1
% gx1=eval(subs(gx,[x y z],{X1, Y1, Z1}))
gx2=subs(gx,[x y z],{X1, Y1, Z1});
f1=subs(f,[x y z],{X1, Y1, Z1});
figure('name','dF/dL vs L');plot(X1,gx2)
figure('name','F vs L');plot(X1,f1)


Y2=2.6,
Z2=0.5:0.05:0.85;
X2=100;
gx2=subs(gx,[x y z],{X2, Y2, Z2})
f2=subs(f,[x y z],{X2, Y2, Z2});
figure('name','dF/d(pitch) vs pitch');plot(Z2,gx2)
figure('name','F vs pitch');plot(Z2,f2)

Y3=1.8:0.1:3.4,
Z3=0.6;
X3=100;
Neff_cap= a/Z3
gx2=subs(gx,[x y z],{X3, Y3, Z3})
f3=subs(f,[x y z],{X3, Y3, Z3});
figure('name','dF/d(Neff) vs Neff');plot(Y3,gx2)
figure('name','F vs Neff');plot(Y3,f3)

Y4=2.6,
Z4=0.6
X4=10:5:100;
Neff_cap= a/Z4
% gx1=eval(subs(gx,[x y z],{X1, Y1, Z1}))
gx2=subs(gx,[x y z],{X4, Y4, Z4});
f4=subs(f,[x y z],{X4, Y4, Z4});
figure('name','dF/dL vs L');plot(X4,gx2)
figure('name','F vs L');plot(X4,f4)

