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
% lam =1.55
% c = 300
% df = 10
% %df = c*(dlam/lam^2)
% dlam = df*(lam^2/c)
%%%%%%%% Start of Copyright Notice %%%%%%%
%
% copyright © 2018
% Written by Sanam Moslemi-Tabrizi
% All rights reserved
% Can NOT be copied or distrutated without including this copyright notice
%
%%%%%%%% End of Copyright Notice %%%%%%%%%%%


% %( c\lambda^2) |delta_lambda| = |delta_f|
% dlam_micron=dlam*(10^-6)
% clear all
% close all

lam =1.55 % micron
d = 0.500 % pitch in micron
dc = 0.6 %duty cycle
n_cl = 1.44
eps_cl=sqrt(n_cl)
eps_s = 1.44^2 % substrate (must be isolator) permittivity
% abs(Neff-lam/d)<1 (1=n_cl?)
Neff_max_1=lam/d + n_cl % Neff < lam/d + n_cl
Neff_min=lam/d - n_cl % Neff > lam/d - n_cl
Neff_max_2 = (2*lam)/d - sqrt(eps_s)
Neff_max = max(Neff_max_1,Neff_max_2)

eps_r=Neff^2
cond1 = (Neff ~= lam/d)

% Neff > sqrt(eps_s)
cond2 = (Neff > sqrt(eps_s))
leakage_factor =((eps_r - eps_cl)*sin(pi*dc))^2;
cond3 = (Neff > n_cl) % close to n_cl


% 
% % lam=1.55
% % nc=1.44   % Top cladding in sio2
% % no=3.4547 % core is Si
% % ne=2.216 % this came from fde with etch-depth = 110 nm
% % z=0
% % r=0.025
% % theta_air=14.5 % in degree
% % % sin(theta_air) = nc * sin(theta_c)
% % theta_c = asind(sind(theta_air)/nc)
% % le = 0.100 % came from design at gcTest.m
% % lo = 0.420 % came from design at gcTest.m
% % pitch = le + lo
% % neff = lam/pitch + sind(theta_air)
% % f0_1=(neff-ne)/(no-ne)
% % f0=lo/pitch
% % f = f0-r*z
% % le=(1-f)*pitch
% % lo=f*pitch
% % neff_1 = f*no+(1-f)*ne
% % theta_air =asind((pitch*neff_1 - lam)/(pitch))
% % thera_c = asind(sind(theta_air)/nc)
% 
% %from the paper: Design Paramet oPTIMZATION OF A SILICON BASED GRATING
% %WAVEGUIDE FOR PREFORMANCE IMPROVEMENT in biochemical sensor application
% %
% %equation (1)
% % the date is chosen by myself but the definations are by the paper
% lam_r=1.55% resonance wavelenght
% ng= 1.60 % group index of waveguide core (core = layer + gratings)
% L= 98 % grating lenght %in the paper the value is fixed at 98 micron
% fwhm_micron= (lam_r^2)/(2*pi*L*ng) % the lenghts are in micron
% fwhm_nano= 1000*((lam_r^2)/(2*pi*L*ng)) % lam_r and L are in micron but fwhm is in nano
% %
% fwhm2= 0.505
% l2=sqrt(2*pi*L*ng*fwhm2/1000)
% %
% fwhm3=0.505
% l3=1.55
% ng3=1000*((l3^2)/(2*pi*L*fwhm3))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % lam = 1.55
% % % % neff= 3.47
% % % % ncl = 1.45
% % % % % n_periods_0 =14
% % % % l_target = 15
% % % % % pitch_0=lam/neff
% % % % pitch = 0.66
% % % % duty_cycle = 0.63
% % % % etch_depth_ends = 0.060
% % % % delta_depth = 0.010
% % % % etch_depth= 0.10
% % % % h_total = 0.22
% % % % etch_mid=h_total - etch_depth
% % % % %
% % % % fill_width=pitch*duty_cycle
% % % % etch_width=pitch*(1-duty_cycle)
% % % % %
% % % % n_first = (etch_depth - etch_depth_ends)/delta_depth
% % % % n_periods_0 = 2*n_first
% % % % L0=n_periods_0 *pitch
% % % % 
% % % % % n_mid = ceil((l_target - L0)/pitch)  % 8
% % % % n_mid=10
% % % % n_extra=2*n_first+n_mid
% % % % 
% % % % n_periods = n_periods_0+n_extra
% % % % l=n_periods *pitch
% % % % 
% % % % for i=1:n_periods
% % % %     x_min_si(i)= pitch*(i-1)+etch_width;
% % % %     x_max_si(i)= pitch*(i);
% % % %     y_min_si(i)= etch_mid;
% % % %     y_max_si(i)= h_total;
% % % % end
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % n_id_1=0;
% % % % pair_id_1=0;
% % % % n_id_2=0;
% % % % pair_id_2=0;
% % % % for i=1:n_periods
% % % %     x_min_etch(i)= pitch*(i);
% % % %     x_max_etch(i)= pitch*(i-1)+etch_width;
% % % %     y_min_si(i)= etch_mid;
% % % %     
% % % %     if (i<=2*n_first)
% % % %         if (n_id_1 >=2)
% % % %             n_id_1=0;
% % % %             pair_id_1=pair_id_1+1;
% % % %         end
% % % %         i
% % % %         etch_depth_n = etch_depth_ends + pair_id_1*delta_depth
% % % %         n_id_1=n_id_1+1;
% % % %     elseif ((i > 2*n_first) && (i <= (2*n_first + n_mid)))
% % % %         i
% % % %         etch_depth_n = etch_depth
% % % %     elseif (i > (2*n_first + n_mid))
% % % %         if (n_id_2 >=2)
% % % %             n_id_2=0;
% % % %             pair_id_2=pair_id_2+1;
% % % %         end
% % % %         i
% % % %          etch_depth_n = etch_depth- (1+pair_id_2)*delta_depth
% % % %          n_id_2=n_id_2+1;
% % % %     end
% % % % end
% % % % 
% % % % 
% % % % % n_periods = ceil(l/pitch)
% % % % 
% % % % %0.03, 0.04 0.05 0.06 0.07 0.08 0.09 0.08 0.07 0.06 0.05 0.04 0.03
% % % % %
% % % % 
% % % % % 1     2     3    4   5    6     7    8   9    10   11    12   13
% % % % 
% % % % 
% % % % theta = asind((pitch*neff-lam)/ncl*pitch)
% % % % % sin(theta) = (pitch*neff-lam)/ncl*pitch