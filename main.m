clear all; close all; clc;

load aeroegm2008.mat

r0 = [2.865408456918535   5.191131097020245   2.848416875743876]*1.0e+06; %m
v0 =  [-5.386247766065933  -0.386715190539288   6.123151881231440]*1.0e+03; %m/s
Tf = 6.218728117616871e+03; %s
tspan = [0 Tf];
N = 120;
M = 120;

tic
[gx, gy, gz] = gravitysphericalharmonic(r0, 'EGM2008',N,'Error');
toc

t_Aero_Toolbox = toc;


tic
%Canonical Unit System
DU = Re; %m
TU = sqrt(DU^3/GM); %s
Req = 1; %Equatorial Radius
mu = 1; %Earth's gravitational parameter

r = r0/DU;
[dRdr,dRdphi,dRdlamda] = dRdr_dRdphi_dRdlamda(N,M,Req,r,mu,C,S);
R = norm(r0);
g2 = -GM/R^3*r0 + (dRdr + dRdphi + dRdlamda)*DU/TU^2;
toc

t_mycode = toc;


Difference = [gx, gy, gz]  - g2
times = t_Aero_Toolbox/t_mycode;

dsiplay = ['The new code ',num2str(times),' times faster'];
disp(dsiplay)

