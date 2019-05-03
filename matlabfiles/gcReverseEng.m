%%%%%%%% Start of Copyright Notice %%%%%%%
%
% copyright © 2018
% Written by Sanam Moslemi-Tabrizi
% All rights reserved
% Can NOT be copied or distrutated without including this copyright notice
%
%%%%%%%% End of Copyright Notice %%%%%%%%%%%
%%  Take from nature's paper: High-efficiency grating-couplers: demonstration of a new design strategy
%% Equations and Data from the paper
% * neff = F*nO + (1-f)*nE
% * pitch = lam  /(neff - sin(theta_air))
% * nc * sin (theta_c)  = sin(theta_air)
% * F = F0 - R*z
% * neff: effective index of the radiative unit
% * nO: effective index for original silicon slab (non etched area)
% * nE: effective index trench (etch area)  % should come from FDE 
% * nC: effective index of top cladding
% * pitch is the period of the grating  = lO + lE
% * lO: length of the not etched area (teeth)
% * lE: length of the etched area (trench)
% * F : fill factor =LO/pitch
% * F0: initial fill factor
% * lam: the coupling wavelength ($\lambda_c$)
% * theta_air: diffraction angle in air
% * theta_c: diffraction angle in top cladding area
% * R: apodization factor
% * z: the distance of each radiative unit from the starting point of the grating
% * h: thickness of the silicon lab (220nm or 260 nm)
% * e: depth of etching
% * lam = $\lambda_c$
%%
%% known variables (from the paper)
%%

lam=1.55 % micron 
nC=1.44  % sio2 in cladding
nO=3.0344 % silicon slab
theta_air=14.5 % in degree
R=0.025

% sin(theta-{air}) = nC * sin(theta_c)
%theta_c = asind(sind(theta_{air})/nC)

%% First radiative unit 
%%
z(1)=0; % distance from the start of the grating to the first radiative unit
lE(1)=0.060; % from table 3 for the paper
lO(1)=0.54;  % from table 3 for the paper
pitch(1)=lE(1) + lO(1);  %from definition of the pitch = fillWidth + echtedWidth
f(1)=lO(1)/pitch(1);  %from definition of the fill factor
%% Bragg condition
% $n_{eff} = \frac{lam}{pitch(1)}+sind(theta_{air})$
%
% * which is equavalent to 
% $n_{eff} = \frac{lam}{pitch(1)}+nC * sind(theta_c)$
%%
%
neff(1)=lam/pitch(1)+sind(theta_air);  % bragg condition for grating
nE=(neff(1)- f(1)*nO)/(1-f(1)) % from equation 1 of the paper. 
%nE is the effective index trench (etch area) but with this set of data we get to a negative value!!
%% Second radiative unit 
%%
z(2)=pitch(1); % distance from the start of the grating to the second radiative unit
lE(2)=0.069; % from table 3 for the paper
lO(2)=0.533; % from table 3 for the paper
pitch(2)=lE(2) + lO(2); %from definition of the pitch = fill-width + echted-width
f(2)=lO(2)/pitch(2); %from definition of the fill factor
f_theo(2)= f(1)-R*z(2); % from equation 3 in the paper
neff(2)=lam/pitch(2)+sind(theta_air);  % bragg condition for grating
nE=(neff(2)- f(2)*nO)/(1-f(2)) % from equation 1 of the paper. 
%nE is the effective index trench (etch area) but with this set of data we get to a negative value!!

