%%%%%%%% Start of Copyright Notice %%%%%%%
%
% copyright © 2018
% Written by Sanam Moslemi-Tabrizi
% All rights reserved
% Can NOT be copied or distrutated without including this copyright notice
%
%%%%%%%% End of Copyright Notice %%%%%%%%%%%
clear all
lam = 1.55
ncore= 2.96
ncl = 1.45
l_target = 15
%pitch = 0.66
duty_cycle = 0.63
%f=lam/pitch;
n_teeth=2.83989
n_groove=2.33587
e=0.12
h=0.22
w=1

neff_theo = duty_cycle*n_teeth+(1-duty_cycle)*n_groove
neff_swg_1_theo = sqrt(duty_cycle*(n_teeth^2)+(1-duty_cycle)*(n_groove^2))
neff_swg_w_theo = 1/(sqrt(duty_cycle*(n_teeth^(-2))+(1-duty_cycle)*(n_groove^(-2))))

% % material_fill = w*h*duty_cycle*pitch*ncore;
% % material_etch= w*(h-e)*(1-duty_cycle)*pitch*ncore+w*(e)*(1-duty_cycle)*pitch*ncl;
% % material_etch_swg_1= w*(1-duty_cycle)*pitch*h*sqrt(((h-e)/h)*ncore^2+(e/h)*ncl^2);
% % material_etch_swg_2= w*(1-duty_cycle)*pitch*h*sqrt(((h-e)/h)*ncore^(-2)+(e/h)*ncl^(-2));
% % v_whole = w*h*pitch;
% % neff_fill_theo = (material_fill)/((duty_cycle)*v_whole)
% % neff_etch_theo = (material_etch)/((1-duty_cycle)*v_whole)
% % neff_etch_theo_swg_1 = (material_etch_swg_1)/((1-duty_cycle)*v_whole)
% % neff_etch_theo_swg_2 = (material_etch_swg_2)/((1-duty_cycle)*v_whole)
% % 
% % 
% % neff_theo = duty_cycle*neff_fill_theo+(1-duty_cycle)*neff_etch_theo
% % neff_theo_swg_1 = sqrt(duty_cycle*neff_fill_theo^2+(1-duty_cycle)*neff_etch_theo^2)
% % neff_theo_swg_2 = 1/sqrt(duty_cycle*neff_fill_theo^(-2)+(1-duty_cycle)*neff_etch_theo^(-2))
% % material_fill_1 = w*h*duty_cycle*pitch*ncore^2;
% % material_etch_1= w*(h-e)*(1-duty_cycle)*pitch*ncore^2+w*(e)*(1-duty_cycle)*pitch*ncl^2;
% % neff_fill_theo_1 = (material_fill_1)/v_whole
% % neff_etch_theo_1 = (material_etch_1)/v_whole
% % neff_theo_1 = sqrt((material_fill_1+material_etch_1)/v_whole)
% % 
% % 
% % material_fill_2 = w*h*duty_cycle*pitch*ncore^(-2);
% % material_etch_2= w*(h-e)*(1-duty_cycle)*pitch*ncore^(-2)+w*(e)*(1-duty_cycle)*pitch*ncl^(-2);
% % neff_fill_theo_2 = (material_fill_2)/v_whole
% % neff_etch_theo_2 = (material_etch_2)/v_whole
% % neff_theo_2 = 1/sqrt((material_fill_2+material_etch_2)/v_whole)
% 
% n_fill_theo = ncore
% n_etch_theo = (e*ncl+(h-e)*ncore)/h
% 
% % etch_depth_ends = 0.060
% % delta_depth = 0.010
% % etch_depth= 0.10
% % h_total = 0.22
% 
% neff=duty_cycle * n_fill + (1-duty_cycle) * n_etch
% neff_swg_parallel = sqrt(duty_cycle * n_fill^2 + (1-duty_cycle) * n_etch^2)
% neff_swg_tang = 1/sqrt(duty_cycle * n_fill^(-2) + (1-duty_cycle) * n_etch^(-2))