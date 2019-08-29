%width_modeSolution
% neff
% material index
close all;
clear all;

materail_neff= cell(50,50);
materail_neff_mod= cell(50,50);


% Read simulation data of 3d fdtd swg_h run.
% data format:  dc_swg - neff mode monitor - neff - width - pitch_swg
swg_h_data_w0p5=load("fdtd_mode_dc_swg_h_multi_w0p5_Data.txt");
swg_h_data_w0p6=load("fdtd_mode_dc_swg_h_multi_w0p6_Data.txt");
swg_h_data_w0p7=load("fdtd_mode_dc_swg_h_multi_w0p7_Data.txt");
swg_h_data_w1p0=load("fdtd_mode_dc_swg_h_multi_w1p0_Data.txt");

% Keep the propagating values
n_swg_h_w0p5=9; % equvalent to dc=0.7
%n_swg_h_w0p5_end=24; % equvalent to dc=0.85
swg_h_data_w0p5=swg_h_data_w0p5(n_swg_h_w0p5:end,:);
n_swg_h_w0p6=8; % equvalent to dc=0.65
swg_h_data_w0p6=swg_h_data_w0p6(n_swg_h_w0p6:end,:);
n_swg_h_w0p7=7; % equvalent to dc=0.6
swg_h_data_w0p7=swg_h_data_w0p7(n_swg_h_w0p7:end,:);


% extract swg_h data
dc_swg_h_w0p5=swg_h_data_w0p5(:,1);
dc_swg_h_w0p6=swg_h_data_w0p6(:,1);
dc_swg_h_w0p7=swg_h_data_w0p7(:,1);
dc_swg_h_w1p0=swg_h_data_w1p0(:,1);
neff_modeMonitor_swg_h_w0p5=swg_h_data_w0p5(:,2);
neff_modeMonitor_swg_h_w0p6=swg_h_data_w0p6(:,2);
neff_modeMonitor_swg_h_w0p7=swg_h_data_w0p7(:,2);
neff_modeMonitor_swg_h_w1p0=swg_h_data_w1p0(:,2);
neff_swg_h_w0p5=swg_h_data_w0p5(:,3);
neff_swg_h_w0p6=swg_h_data_w0p6(:,3);
neff_swg_h_w0p7=swg_h_data_w0p7(:,3);
neff_swg_h_w1p0=swg_h_data_w1p0(:,3);
width_swg_h_w0p5=swg_h_data_w0p5(1,4);
width_swg_h_w0p6=swg_h_data_w0p6(1,4);
width_swg_h_w0p7=swg_h_data_w0p7(1,4);
width_swg_h_w1p0=swg_h_data_w1p0(1,4);


% Read simulation data of 3d fdtd swg_v run.
% data format:  dc_swg - neff - width - pitch_swg
swg_v_data_w0p5=load("fdtd_mode_dc_swg_multi_w0p5_Data.txt");
swg_v_data_w0p6=load("fdtd_mode_dc_swg_multi_w0p6_Data.txt");
%swg_v_data_w0p7=load("fdtd_mode_dc_swg_multi_w0p5_Data.txt");
% Keep the propagating values
n_swg_v_w0p5=11; % equvalent to dc=0.5
swg_v_data_w0p5=swg_v_data_w0p5(1:n_swg_v_w0p5,:);
n_swg_v_w0p6=13; % equvalent to dc=0.4
swg_v_data_w0p6=swg_v_data_w0p6(1:n_swg_v_w0p6,:);



% extract swg_v data
dc_swg_v_w0p5=swg_v_data_w0p5(:,1);
dc_swg_v_w0p6=swg_v_data_w0p6(:,1);
%dc_swg_v_w0p7=swg_v_data_w0p7(:,1);
neff_swg_v_w0p5=swg_v_data_w0p5(:,2);
neff_swg_v_w0p6=swg_v_data_w0p6(:,2);
%neff_swg_v_w0p7=swg_v_data_w0p7(:,2);
width_swg_v_w0p5=swg_v_data_w0p5(1,3);
width_swg_v_w0p6=swg_v_data_w0p6(1,3);
%width_swg_v_w0p7=swg_v_data_w0p7(1,3);


% Read simulation data of 3d fdtd for each effMaterial run.
% data format:  index of the under test effMaterial -  neff mode monitor -  neff 3d fdtd -  width 
fdtddata=load("fdtd_effMaterial_Index_effMaterial_multi_Data.txt");
fdtddata1=load("fdtd_effMaterial_Index_effMaterial1_multi_Data.txt");
fdtddata2=load("fdtd_effMaterial_Index_effMaterial2_multi_Data.txt");
fdtddata3=load("fdtd_effMaterial_Index_effMaterial3_multi_Data_MOD.txt");
fdtddata4=load("fdtd_effMaterial_Index_effMaterial4_multi_Data_MOD.txt");
fdtddata5=load("fdtd_effMaterial_Index_effMaterial5_multi_Data_MOD.txt");
% extract effMaterial 3d FDTD data
neff_3d_fdtd=fdtddata(:,3);
neff_3d_fdtd1=fdtddata1(:,3);
neff_3d_fdtd2=fdtddata2(:,3);
neff_3d_fdtd3=fdtddata3(:,3);
neff_3d_fdtd4=fdtddata4(:,3);
neff_3d_fdtd5=fdtddata5(:,3);
%
materialIndex=fdtddata(:,1);  % = effMaterialIndex*ones(n,1)
materialIndex1=fdtddata1(:,1); % = effMaterialIndex1*ones(n,1)
materialIndex2=fdtddata2(:,1); % = effMaterialIndex2*ones(n,1)
materialIndex3=fdtddata3(:,1); % = effMaterialIndex3*ones(n,1)
materialIndex4=fdtddata4(:,1); % = effMaterialIndex4*ones(n,1)
materialIndex5=fdtddata5(:,1); % = effMaterialIndex5*ones(n,1)
%
width_3d_fdtd=fdtddata(:,4); % = waveguideWidth*ones(n,1)
width_3d_fdtd1=fdtddata1(:,4); % = waveguideWidth*ones(n,1)
width_3d_fdtd2=fdtddata2(:,4); % = waveguideWidth*ones(n,1)
width_3d_fdtd3=fdtddata3(:,4); % = waveguideWidth*ones(n,1)
width_3d_fdtd4=fdtddata4(:,4); % = waveguideWidth*ones(n,1)
width_3d_fdtd5=fdtddata5(:,4); % = waveguideWidth*ones(n,1)
% NOTE width_3d_fdtd = width_3d_fdtd1= width_3d_fdtd2= width_3d_fdtd3 = width_3d_fdtd4 = width_3d_fdtd5
%
tf_width = isequal(width_3d_fdtd,width_3d_fdtd1,width_3d_fdtd2,width_3d_fdtd3,width_3d_fdtd4,width_3d_fdtd5);
if (tf_width==0)
    disp('ERROR: with WIDTH sizes')
end

% extra neff_Materail 3D FDFD simulation data for w0p7
% 1:material index , 2:neff (mode)
extra_mat_neff_w0p7=load('extra3_fdtd_effMaterial_Index_multi_w0p7_Data.txt');


% extra neff_Materail 3D FDFD simulation data for w0p6 
extra_mat_w0p6=[2.5342; 2.6224; 2.6866; 2.7738; 2.8702 ];
extra_mat_neff_w0p6=[2.5342 , 1.58445 ; 2.6224,1.65227 ; 2.6866,1.70555 ; 2.7738,1.7851 ; 2.8702,1.8808 ];
%

%extra_dc_swg_h = 0.73, 0.78, 0.82, 0.87 3D FDFD simulation data for w0p5 
extra_dc_swg_h_w0p5_Data=load('fdtd_mode_dc_swg_h_multi_w0p5_Data_EXTRA.txt'); %1=dc, 2=mode monitor, 3=neff, 4=width, 5=pitch
%

%extra_dc_swg_h = 0.68, 0.73, 0.78, 0.82 3D FDFD simulation data for w0p6 
extra_dc_swg_h_w0p6_Data=load('fdtd_mode_dc_swg_h_multi_w0p6_Data_EXTRA.txt');%1=dc, 2=mode monitor, 3=neff, 4=width, 5=pitch 
%

%extra_dc_swg_h = 0.68, 0.73, 0.78, 0.82 3D FDFD simulation data for w0p7 
extra_dc_swg_h_w0p7_Data=load('fdtd_mode_dc_swg_h_multi_w0p7_Data_EXTRA.txt');%1=dc, 2=mode monitor, 3=neff, 4=width, 5=pitch 
%

%extra_dc_swg_h_H =  0.757 , 0.802,  0.833, 0.874,  0.915. the values of dc
%are the fitted values for the extra_mat_neff_w0p6. We need to put these
%set of datapoints on the fitted_mode_material_vs_dc_swg_h_w0p6 graph CASE A
extra_dc_swg_h_H_w0p6_Data=load('fitted_dc_swg_h_H_Data.txt'); 
%

%extra_dc_swg_h_H_w0p7 =  0.645, 0.691, 0.748, 0.793, 0.829 0.868 the values of dc
%are the fitted values for the extra_mat_neff_w0p7. We need to put these
%set of datapoints on the fitted_mode_material_vs_dc_swg_h_w0p7 graph CASE A
%extra_dc_swg_h_H_w0p7_Data=load('fitted_dc_swg_h_H_w0p7_Data.txt'); 
extra_dc_swg_h_H_w0p7_Data=load('fitted2_dc_swg_h_H_w0p7_Data.txt'); 
%

%extra_dc_swg_v_H =  0.5507, 0.5982, 0.6273, 0.6661, 0.7156. the values of dc
%are the fitted values for the extra_mat_neff_w0p6. We need to put these
%set of datapoints on the itted_mode_material_vs_dc_swg_v_w0p6 graph CASE B
extra_dc_swg_v_H_w0p6_Data=load('fitted_dc_swg_v_H_Data.txt'); 


% creat a matrix of 3d fdtd data
% dataformat wgWidth - neff for effMaterial - neff for effMaterial1 - neff for effMaterial2 - neff for effMaterial3 - neff for effMaterial4 - neff for effMaterial5
% for lower index materials in lower width sizes we don't have propagaion
% and I put neff=zero for those elements.
mat_neff_3d_fdtd=[width_3d_fdtd,neff_3d_fdtd,neff_3d_fdtd1,neff_3d_fdtd2,neff_3d_fdtd3,neff_3d_fdtd4,neff_3d_fdtd5];
%
% an array of material indecis 3.47 , 3.25 , 3 , 2.75, 2.5 , 2.25
material_index_range=[materialIndex(1,1), materialIndex1(1,1), materialIndex2(1,1), materialIndex3(1,1), materialIndex4(1,1), materialIndex5(1,1)];


% remove the data for lower index material 92.25) and low wg width (0.35, 0.4)
mat_neff_3d_fdtd_mod=mat_neff_3d_fdtd(3:end,1:6);
material_index_range_mod=material_index_range(1,1:5);
w_mod=width_3d_fdtd(3:end);
% Modifed wg width: 0.45, 0.5,0.55, ..1.6, 1.65, 1.7


% for each wg_width (i) create a matrix of material index and the resulted mode index (neff) from the 3d fdtd simulations
% 3d FDFD sims on effMaterial
% data format for each width: MaterialIndex - ModeIndex(neff)
for i=1:length(mat_neff_3d_fdtd_mod) 
    materail_neff_mod{i}=[material_index_range_mod',mat_neff_3d_fdtd_mod(i,2:end)'];
end
material_range=materail_neff_mod{i}(:,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CURVE FITTING %%%%%%%%%%%%%
% ==============> SWG_H curve fitting
% for dc_swg_h=0.7 ==> start propagating , neff > n_cladding , SWG_H
%n_swg_h_w0p5=9;
fit_mode_vs_dc_swg_h_w0p5=fit(swg_h_data_w0p5(:,1),swg_h_data_w0p5(:,3),'poly4'); % To get fitted neff for given dc swg_h
%fit_dc_swg_h_vs_mode_w0p5=fit(swg_h_data_w0p5(:,3),swg_h_data_w0p5(:,1),'poly4'); % To get fitted dc swg_h for given neff
xdata_swg_h_w0p5=swg_h_data_w0p5(:,3);
ydata_swg_h_w0p5=swg_h_data_w0p5(:,1);

fit_dc_swg_h_vs_mode_w0p5 = @(b,xdata_swg_h_w0p5) b(1).*exp(b(2).*xdata_swg_h_w0p5)+b(3);  % Objective Function
B1 = fminsearch(@(b) norm(ydata_swg_h_w0p5 - fit_dc_swg_h_vs_mode_w0p5(b,xdata_swg_h_w0p5)), [-200; -1; 100]);  
%
% for dc_swg_h=0.65 ==> start propagating , neff > n_cladding , SWG_H
%n_swg_h_w0p6=8;
fit_mode_vs_dc_swg_h_w0p6=fit(swg_h_data_w0p6(:,1),swg_h_data_w0p6(:,3),'poly4'); % To get fitted neff for given dc swg_h
%fit_dc_swg_h_vs_mode_w0p6=fit(swg_h_data_w0p6(:,3),swg_h_data_w0p6(:,1),'poly4'); % To get fitted dc swg_h for given neff
xdata_swg_h_w0p6=swg_h_data_w0p6(:,3);
ydata_swg_h_w0p6=swg_h_data_w0p6(:,1);

fit_dc_swg_h_vs_mode_w0p6 = @(b,xdata_swg_h_w0p6) b(1).*exp(b(2).*xdata_swg_h_w0p6)+b(3);  % Objective Function
B2 = fminsearch(@(b) norm(ydata_swg_h_w0p6 - fit_dc_swg_h_vs_mode_w0p6(b,xdata_swg_h_w0p6)), [-200; -1; 100]);  
%
% for dc_swg_h=0.6 ==> start propagating , neff > n_cladding , SWG_H
%n_swg_h_w0p7=7;
fit_mode_vs_dc_swg_h_w0p7=fit(swg_h_data_w0p7(:,1),swg_h_data_w0p7(:,2),'poly4'); % To get fitted neff for given dc swg_h
%fit_dc_swg_h_vs_mode_w0p7=fit(swg_h_data_w0p7(:,2),swg_h_data_w0p7(:,1),'poly6'); % To get fitted dc swg_h for given neff
xdata_swg_h_w0p7=swg_h_data_w0p7(:,3);
ydata_swg_h_w0p7=swg_h_data_w0p7(:,1);

fit_dc_swg_h_vs_mode_w0p7 = @(b,xdata_swg_h_w0p7) b(1).*exp(b(2).*xdata_swg_h_w0p7)+b(3);  % Objective Function
B3 = fminsearch(@(b) norm(ydata_swg_h_w0p7 - fit_dc_swg_h_vs_mode_w0p7(b,xdata_swg_h_w0p7)), [-200; -1; 100]); 
%
% ==============> SWG_V curve fitting
% for w=0.5 ==> id=2 SWG_V
fit_mode_vs_dc_swg_v_w0p5=fit(swg_v_data_w0p5(:,1),swg_v_data_w0p5(:,2),'poly6'); % To get fitted neff for given dc swg_v
fit_dc_swg_v_vs_mode_w0p5=fit(swg_v_data_w0p5(:,2),swg_v_data_w0p5(:,1),'poly6'); % To get fitted dc swg_v for given neff
%
% for w=0.6 ==> id=4 SWG_V
fit_mode_vs_dc_swg_v_w0p6=fit(swg_v_data_w0p6(:,1),swg_v_data_w0p6(:,2),'poly6'); % To get fitted neff for given dc swg_v
fit_dc_swg_v_vs_mode_w0p6=fit(swg_v_data_w0p6(:,2),swg_v_data_w0p6(:,1),'poly6'); % To get fitted dc swg_v for given neff
%
% ==============>  3D FDTD curve fitting
% for w=0.5 ==> id=2 % simulated data
w_f_w0p5=w_mod(2);
f_w0p5=materail_neff_mod{2};  
fit_mode_vs_material_w0p5=fit(f_w0p5(:,1),f_w0p5(:,2),'poly4'); % To get fitted neff for given effMaterial Index
fit_material_vs_mode_w0p5=fit(f_w0p5(:,2),f_w0p5(:,1),'poly4'); % To get fitted effMaterial Index for given neff 
%
% for w=0.6 ==> id=4
w_f_w0p6=w_mod(4);
f_w0p6=materail_neff_mod{4}; % simulated data
fit_mode_vs_material_w0p6=fit(f_w0p6(:,1),f_w0p6(:,2),'poly4'); % To get fitted neff for given effMaterial Index
fit_material_vs_mode_w0p6=fit(f_w0p6(:,2),f_w0p6(:,1),'poly4'); % To get fitted effMaterial Index for given neff
%
% for w=0.7 ==> id=6
w_f_w0p7=w_mod(6);
f_w0p7=materail_neff_mod{6}; % simulated data
fit_mode_vs_material_w0p7=fit(f_w0p7(:,1),f_w0p7(:,2),'poly4'); % To get fitted neff for given effMaterial Index
fit_material_vs_mode_w0p7=fit(f_w0p7(:,2),f_w0p7(:,1),'poly4'); % To get fitted effMaterial Index for given neff 
%

% % % % % % % %  REDUNDANT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REDUNDANT FULL CURVE FITTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % 
% % % % % % % % FULL CURVE FITTING swg_v
% % % % % % % % w0p5
% % % % % % % % Input is dc_swg_v => Calculate the fitted neff => Calculate the fitted effMaterialIndex
% % % % % % % neff_fitted_w0p5_1=fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5);
% % % % % % % material_fitted_w0p5_1=fit_material_vs_mode_w0p5(neff_fitted_w0p5_1);
% % % % % % % %
% % % % % % % % w0p6
% % % % % % % % Input is dc_swg_v => Calculate the fitted neff => Calculate the fitted effMaterialIndex
% % % % % % % neff_fitted_w0p6_1=fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6);
% % % % % % % material_fitted_w0p6_1=fit_material_vs_mode_w0p6(neff_fitted_w0p6_1);
% % % % % % % %
% % % % % % % 
% % % % % % % 
% % % % % % % 
% % % % % % % %  REDUNDANT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   REDUNDANT EXTRA DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % %%%%% EXTRA DATA for SWG_V W0p6  %%%%%
% % % % % % % extra_fitted_dc_swg_v_vs_neff_w0p6=fit_dc_swg_v_vs_mode_w0p6(extra_mat_neff_w0p6(:,2))
% % % % % % % extra_material_vs_neff_w0p6=fit_material_vs_mode_w0p6(extra_mat_neff_w0p6(:,2))
% % % % % % % extra_fitted_neff_vs_material_w0p6=fit_mode_vs_material_w0p6(extra_material_vs_neff_w0p6)
% % % % % % % %
% % % % % % % %%%%% EXTRA DATA for SWG_H W0p6  %%%%%
% % % % % % % extra_fitted_dc_swg_h_vs_neff_w0p6=fit_dc_swg_h_vs_mode_w0p6(B2,extra_mat_neff_w0p6(:,2))
% % % % % % % extra_material_vs_neff_w0p6_2=fit_material_vs_mode_w0p6(extra_mat_neff_w0p6(:,2))
% % % % % % % extra_fitted_neff_vs_material_w0p6_2=fit_mode_vs_material_w0p6(extra_material_vs_neff_w0p6)
% % % % % % % %

%%%  Scenario 1 %%%% Full process for w0p6
%Input is material_H (extra_mat_neff_w0p6) => sim_neff_h , fitted_neff_h
% Put is on material_and_mode_vs_dc_swg  plots for w0p6 
%Output is dc_swg_v_H
%Output is dc_swg_h_H
material_H=extra_mat_neff_w0p6(:,1);
sim_neff_H=extra_mat_neff_w0p6(:,2);
fitted_neff_H=fit_mode_vs_material_w0p6(material_H);
error_neff_H_material=100*(sim_neff_H-fitted_neff_H)./sim_neff_H;
fitted_dc_swg_v_H_vs_mode_w0p6=fit_dc_swg_v_vs_mode_w0p6(sim_neff_H);
% % fitted_dc_swg_h_H_vs_mode_w0p6=fit_dc_swg_h_vs_mode_w0p6(B2,sim_neff_H);
%%%%%%% NOTE extra_dc_swg_v_H_w0p6_Data(:,1) == fitted_dc_swg_v_H_vs_mode_w0p6
% % %%%%%%% NOTE extra_dc_swg_h_H_w0p6_Data(:,1) == fitted_dc_swg_h_H_vs_mode_w0p6
% % sim_neff_H_swg_h= extra_dc_swg_h_H_w0p6_Data(:,3);
sim_neff_H_swg_v= extra_dc_swg_v_H_w0p6_Data(:,2);
% % fitted_neff_H_swg_h=fit_mode_vs_dc_swg_h_w0p6(extra_dc_swg_h_H_w0p6_Data(:,1));
fitted_neff_H_swg_v=fit_mode_vs_dc_swg_v_w0p6(extra_dc_swg_v_H_w0p6_Data(:,1));
% % fitted_material_H_swg_h=fit_material_vs_mode_w0p6(sim_neff_H_swg_h);
fitted_material_H_swg_v=fit_material_vs_mode_w0p6(sim_neff_H_swg_v);
% % error_fitted_material_H_swg_h=100*(material_H-fitted_material_H_swg_h)./material_H
error_fitted_material_H_swg_v=100*(material_H-fitted_material_H_swg_v)./material_H
%


%%% Scenario 1 %%%% Full process for w0p7 JUST FOR W0P7 SWG_H
%Input is material_H (extra_mat_neff_w0p7) => sim_neff_h , fitted_neff_h
% Put is on  plots for w0p7 (it is already done)
%Output is dc_swg_h_H

material_H_w0p7=extra_mat_neff_w0p7(:,1);
sim_neff_H_w0p7=extra_mat_neff_w0p7(:,3);
fitted_neff_H_w0p7=fit_mode_vs_material_w0p7(material_H_w0p7);
error_neff_H_material_w0p7=100*(sim_neff_H_w0p7-fitted_neff_H_w0p7)./sim_neff_H_w0p7
fitted_dc_swg_h_H_vs_mode_w0p7=fit_dc_swg_h_vs_mode_w0p7(B3,sim_neff_H_w0p7)
%%%%%%% NOTE extra_dc_swg_h_H_w0p7_Data(:,1) == fitted_dc_swg_h_H_vs_mode_w0p7
sim_neff_H_swg_h_w0p7= extra_dc_swg_h_H_w0p7_Data(:,3);
fitted_neff_H_swg_h_w0p7=fit_mode_vs_dc_swg_h_w0p7(extra_dc_swg_h_H_w0p7_Data(:,1));
fitted_material_H_swg_h_w0p7=fit_material_vs_mode_w0p7(sim_neff_H_swg_h_w0p7);
error_fitted_material_H_swg_h_w0p7=100*(material_H_w0p7-fitted_material_H_swg_h_w0p7)./material_H_w0p7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Scenario 2 %%% Full process for w0p6 SWG_V
%Input is dc_swg_v (extra_dc_swg_v_w0p6_Data) => sim_neff_L , fitted_neff_L
% Put is on  plots for w0p6 
%Output is material_h_L
dc_swg_v_L=extra_dc_swg_v_H_w0p6_Data(:,1);
sim_neff_swg_v_L=extra_dc_swg_v_H_w0p6_Data(:,2);
fitted_neff_swg_v_L=fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_L);
fitted_material_swg_v_L_w0p6=fit_material_vs_mode_w0p6(sim_neff_swg_v_L);
%%%%%%% NOTE material_H_w0p6_Data(:,1) == fitted_material_swg_v_L_w0p6
sim_neff_v_L_material_w0p6=fit_mode_vs_material_w0p6(fitted_material_swg_v_L_w0p6)



%%% Scenario 2 %%% Full process for w0p7 SWG_H
%Input is dc_swg_h (extra_dc_swg_h_w0p7_Data) => sim_neff_L , fitted_neff_L
% Put is on  plots for w0p7 
%Output is material_h_L
dc_swg_h_L=extra_dc_swg_h_H_w0p7_Data(:,1);
sim_neff_swg_h_L=extra_dc_swg_h_H_w0p7_Data(:,3);
fitted_neff_swg_h_L=fit_mode_vs_dc_swg_h_w0p7(dc_swg_h_L);
fitted_material_swg_h_L_w0p7=fit_material_vs_mode_w0p7(sim_neff_swg_h_L);
%%%%%%% NOTE material_H_w0p7(:,1) == fitted_material_swg_h_L_w0p7
sim_neff_h_L_material_w0p7=fit_mode_vs_material_w0p7(fitted_material_swg_h_L_w0p7)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('PLOTING SWG_H GRAPHS W0.5')
figure('name','dc_swg_h_vs_mode_w0p5_sim_fitted');
plot(neff_swg_h_w0p5(:),dc_swg_h_w0p5(:));
hold on
plot(swg_h_data_w0p5(:,3),fit_dc_swg_h_vs_mode_w0p5(B1,swg_h_data_w0p5(:,3)));
legend('sim','fitted')
title('dc\_swg\_h\_vs\_mode\_w0p5\_sim\_fitted');
xlabel('Neff(mode index)');
ylabel('swg h duty cycle, sim and fitted');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
grid on


figure('name','mode_vs_dc_swg_h_w0p5_sim_fitted');
plot(dc_swg_h_w0p5(:),neff_swg_h_w0p5(:));
hold on
plot(dc_swg_h_w0p5(:),fit_mode_vs_dc_swg_h_w0p5(dc_swg_h_w0p5(:)))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_h\_w0p5\_sim\_fitted');
xlabel('swg h duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
%

disp('PLOTING SWG_H GRAPHS W0.6')
figure('name','dc_swg_h_vs_mode_w0p6_sim_fitted');
plot(neff_swg_h_w0p6(:),dc_swg_h_w0p6(:));
hold on
plot(swg_h_data_w0p6(:,3),fit_dc_swg_h_vs_mode_w0p6(B2,swg_h_data_w0p6(:,3)));
legend('sim','fitted')
title('dc\_swg\_h\_vs\_mode\_w0p6\_sim\_fitted');
xlabel('Neff(mode index)');
ylabel('swg h duty cycle, sim and fitted');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);


figure('name','mode_vs_dc_swg_h_w0p6_sim_fitted');
plot(dc_swg_h_w0p6(:),neff_swg_h_w0p6(:));
hold on
plot(dc_swg_h_w0p6(:),fit_mode_vs_dc_swg_h_w0p6(dc_swg_h_w0p6(:)))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_h\_w0p6\_sim\_fitted');
xlabel('swg h duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);
%
 
disp('PLOTING SWG_H GRAPHS W0.7')
figure('name','dc_swg_h_vs_mode_w0p7_sim_fitted');
plot(neff_swg_h_w0p7(:),dc_swg_h_w0p7(:));
hold on
plot(swg_h_data_w0p7(:,3),fit_dc_swg_h_vs_mode_w0p7(B3,swg_h_data_w0p7(:,3)));
legend('sim','fitted')
title('dc\_swg\_h\_vs\_mode\_w0p7\_sim\_fitted');
xlabel('Neff(mode index)');
ylabel('swg h duty cycle, sim and fitted');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);


figure('name','mode_vs_dc_swg_h_w0p7_sim_fitted');
plot(dc_swg_h_w0p7(:),neff_swg_h_w0p7(:));
hold on
plot(dc_swg_h_w0p7(:),fit_mode_vs_dc_swg_h_w0p7(dc_swg_h_w0p7(:)))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_h\_w0p7\_sim\_fitted');
xlabel('swg h duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.7 micron'};
text(0.5,3, dims);
%
% =================================================================

disp('PLOTING SWG_V GRAPHS W0.5')
figure('name','dc_swg_v_vs_mode_w0p5_sim_fitted');
plot(neff_swg_v_w0p5,dc_swg_v_w0p5);
hold on
plot(fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5),dc_swg_v_w0p5)
legend('sim','fitted')
title('dc\_swg\_v\_vs\_mode\_w0p5\_sim\_fitted');
xlabel('Neff(mode index)');
ylabel('swg v duty cycle, sim and fitted');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);

figure('name','mode_vs_dc_swg_v_w0p5_sim_fitted');
plot(dc_swg_v_w0p5,neff_swg_v_w0p5);
hold on
plot(dc_swg_v_w0p5,fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_v\_w0p5\_sim\_fitted');
xlabel('swg v duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
%

disp('PLOTING SWG_V GRAPHS W0.6')
figure('name','dc_swg_v_vs_mode_w0p6_sim_fitted');
plot(neff_swg_v_w0p6,dc_swg_v_w0p6);
hold on
plot(fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6),dc_swg_v_w0p6)
legend('sim','fitted')
title('dc\_swg\_v\_vs\_mode\_w0p6\_sim\_fitted');
xlabel('Neff(mode index)');
ylabel('swg v duty cycle, sim and fitted');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);

figure('name','mode_vs_dc_swg_v_w0p6_sim_fitted');
plot(dc_swg_v_w0p6,neff_swg_v_w0p6);
hold on
plot(dc_swg_v_w0p6,fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_v\_w0p6\_sim\_fitted');
xlabel('swg v duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);
% =================================================================

disp('PLOTING EffMaterial GRAPHS W0.5')
figure('name','material_vs_neff_w0p5');
plot(f_w0p5(:,2),f_w0p5(:,1),'r',f_w0p5(:,2),fit_material_vs_mode_w0p5(f_w0p5(:,2)),'*b');
title('MaterialIndex\_ModeIndex\_w0p5\_sims')
xlabel('Mode Index');
ylabel('Material Index');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
legend('sim','fitted');



figure('name','neff_vs_material_w0p5');
plot(f_w0p5(:,1),f_w0p5(:,2),'r',f_w0p5(:,1),fit_mode_vs_material_w0p5(f_w0p5(:,1)),'*b');
title('ModeIndex\_MaterialIndex\_w0p5\_sims');
xlabel('Material Index');
ylabel('Mode Index');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
legend('sim','fitted');
%


disp('PLOTING EffMaterial GRAPHS W0.6')  % IMPORTANT GRAPH
figure('name','SMT_material_vs_neff_w0p6');
plot(f_w0p6(:,2),f_w0p6(:,1),'-.r',f_w0p6(:,2),fit_material_vs_mode_w0p6(f_w0p6(:,2)),'--b');
hold on
plot(extra_mat_neff_w0p6(1,2),extra_mat_neff_w0p6(1,1),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_fitted_neff_vs_material_w0p6(1,1),extra_material_vs_neff_w0p6(1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(2,2),extra_mat_neff_w0p6(2,1),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_fitted_neff_vs_material_w0p6(2,1),extra_material_vs_neff_w0p6(2),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(3,2),extra_mat_neff_w0p6(3,1),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_fitted_neff_vs_material_w0p6(3,1),extra_material_vs_neff_w0p6(3),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(4,2),extra_mat_neff_w0p6(4,1),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_fitted_neff_vs_material_w0p6(4,1),extra_material_vs_neff_w0p6(4),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(5,2),extra_mat_neff_w0p6(5,1),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_fitted_neff_vs_material_w0p6(5,1),extra_material_vs_neff_w0p6(5),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
legend('sim','fitted','extra sim1','extra sim2','extra sim3','extra sim4','extra sim5');
title('MaterialIndex\_ModeIndex\_w0p6\_sims');
xlabel('Mode Index');
ylabel('Material Index');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);

figure('name','SMT_neff_vs_material_w0p6'); %  IMPORTANT GRAPH
plot(f_w0p6(:,1),f_w0p6(:,2),'-.r',f_w0p6(:,1),fit_mode_vs_material_w0p6(f_w0p6(:,1)),'--b');
hold on
plot(extra_mat_neff_w0p6(1,1),extra_mat_neff_w0p6(1,2),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_material_vs_neff_w0p6(1),extra_fitted_neff_vs_material_w0p6(1,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(2,1),extra_mat_neff_w0p6(2,2),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_material_vs_neff_w0p6(2),extra_fitted_neff_vs_material_w0p6(2,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(3,1),extra_mat_neff_w0p6(3,2),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_material_vs_neff_w0p6(3),extra_fitted_neff_vs_material_w0p6(3,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(4,1),extra_mat_neff_w0p6(4,2),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_material_vs_neff_w0p6(4),extra_fitted_neff_vs_material_w0p6(4,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
plot(extra_mat_neff_w0p6(5,1),extra_mat_neff_w0p6(5,2),'-p','MarkerEdgeColor','green','MarkerSize',10,'LineWidth',2)
%plot(extra_material_vs_neff_w0p6(5),extra_fitted_neff_vs_material_w0p6(5,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)

%
legend('sim','fitted','extra sim1','extra sim2','extra sim3','extra sim4','extra sim5');
title('ModeIndex\_MaterialIndex\_w0p6\_sims');
xlabel('Material Index');
ylabel('Mode Index');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);
%

disp('PLOTING EffMaterial GRAPHS W0.7')
figure('name','material_vs_neff_w0p7');
plot(f_w0p7(:,2),f_w0p7(:,1),'r',f_w0p7(:,2),fit_material_vs_mode_w0p7(f_w0p7(:,2)),'*b');
title('MaterialIndex\_ModeIndex\_w0p7\_sims');
xlabel('Mode Index');
ylabel('Material Index');
dims={'wg width=0.7 micron'};
text(0.5,3, dims);
legend('sim','fitted');

figure('name','neff_vs_material_w0p7');
plot(f_w0p7(:,1),f_w0p7(:,2),'r',f_w0p7(:,1),fit_mode_vs_material_w0p7(f_w0p7(:,1)),'*b');
title('ModeIndex\_MaterialIndex\_w0p7\_sims');
xlabel('Material Index');
ylabel('Mode Index');
dims={'wg width=0.7 micron'};
text(0.5,3, dims);
legend('sim','fitted');

% % SWG_V
% disp('PLOTING Mode and Material vs DC SWG_V W0.5')
% figure('name','fitted_mode_material_vs_dc_swg_v_w0p5');
% %'fitted neff','fitted effMaterial for fitted neff'
% plot(dc_swg_v_w0p5,fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5),'r',dc_swg_v_w0p5,fit_material_vs_mode_w0p5(fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5)),'b')
% hold on
% %'sim neff','fitted effmaterial for sim neff'
% plot(dc_swg_v_w0p5,neff_swg_v_w0p5,'*r',dc_swg_v_w0p5,fit_material_vs_mode_w0p5(neff_swg_v_w0p5),'*b')
% title('Mode\_MaterialIndex\_vs\_dc\_swg\_v\_w0p5\_sims');
% xlabel('dc swg v');
% ylabel('Mode Index, Material Index');
% dims={'wg width=0.5 micron'};
% text(0.5,3, dims);
% legend('fitted neff (input:dc\_swg\_v)','fitted effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)')
% %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%
% scenario 1

% Choose Material Index
% find simulated neff
% find fitted neff
% find fitted dc_swg
% run 3d fdtd with the fitted dc_swg and get the simulated neff


% SWG_V
disp('PLOTING Mode and Material vs DC SWG_V W0.6')
figure('name','SMT_fitted_mode_material_vs_dc_swg_v_w0p6');  % IMPORTANT GRAPH
plot(dc_swg_v_w0p6,fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6),'-.r',dc_swg_v_w0p6,fit_material_vs_mode_w0p6(fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6)),'-.b')
hold on
plot(dc_swg_v_w0p6,neff_swg_v_w0p6,'--c',dc_swg_v_w0p6,fit_material_vs_mode_w0p6(neff_swg_v_w0p6),'--m')
plot(extra_dc_swg_v_H_w0p6_Data(1,1),sim_neff_H_swg_v(1,1),'^k',extra_dc_swg_v_H_w0p6_Data(1,1),fitted_material_H_swg_v(1,1),'^k')
plot(extra_dc_swg_v_H_w0p6_Data(2,1),sim_neff_H_swg_v(2,1),'^k',extra_dc_swg_v_H_w0p6_Data(2,1),fitted_material_H_swg_v(2,1),'^k')
plot(extra_dc_swg_v_H_w0p6_Data(3,1),sim_neff_H_swg_v(3,1),'^k',extra_dc_swg_v_H_w0p6_Data(3,1),fitted_material_H_swg_v(3,1),'^k')
plot(extra_dc_swg_v_H_w0p6_Data(4,1),sim_neff_H_swg_v(4,1),'^k',extra_dc_swg_v_H_w0p6_Data(4,1),fitted_material_H_swg_v(4,1),'^k')
plot(extra_dc_swg_v_H_w0p6_Data(5,1),sim_neff_H_swg_v(5,1),'^k',extra_dc_swg_v_H_w0p6_Data(5,1),fitted_material_H_swg_v(5,1),'^k')
plot(fitted_dc_swg_v_H_vs_mode_w0p6(1,1),sim_neff_H(1,1),'^g',fitted_dc_swg_v_H_vs_mode_w0p6(1,1),material_H(1,1),'^g')
plot(fitted_dc_swg_v_H_vs_mode_w0p6(2,1),sim_neff_H(2,1),'^g',fitted_dc_swg_v_H_vs_mode_w0p6(2,1),material_H(2,1),'^g')
plot(fitted_dc_swg_v_H_vs_mode_w0p6(3,1),sim_neff_H(3,1),'^g',fitted_dc_swg_v_H_vs_mode_w0p6(3,1),material_H(3,1),'^g')
plot(fitted_dc_swg_v_H_vs_mode_w0p6(4,1),sim_neff_H(4,1),'^g',fitted_dc_swg_v_H_vs_mode_w0p6(4,1),material_H(4,1),'^g')
plot(fitted_dc_swg_v_H_vs_mode_w0p6(5,1),sim_neff_H(5,1),'^g',fitted_dc_swg_v_H_vs_mode_w0p6(5,1),material_H(5,1),'^g')
title('Mode\_MaterialIndex\_vs\_dc\_swg\_v\_w0p6\_sims');
xlabel('dc swg v');
ylabel('Mode Index, Material Index');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);
legend('fitted neff (input:dc\_swg\_v)','effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)')
%
% % % % SWG_H   ???????????? NEED TO REPEAT FOR W0P7
% % % disp('PLOTING Mode and Material vs DC SWG_H W0.6')
% % % figure('name','fitted_mode_material_vs_dc_swg_h_w0p6');
% % % plot(dc_swg_h_w0p6,fit_mode_vs_dc_swg_h_w0p6(dc_swg_h_w0p6),'r',dc_swg_h_w0p6,fit_material_vs_mode_w0p6(fit_mode_vs_dc_swg_h_w0p6(dc_swg_h_w0p6)),'b')
% % % hold on
% % % plot(dc_swg_h_w0p6,neff_swg_h_w0p6,'*r',dc_swg_h_w0p6,fit_material_vs_mode_w0p6(neff_swg_h_w0p6),'*b')
% % % plot(extra_dc_swg_h_H_w0p6_Data(1,1),sim_neff_H_swg_h(1,1),'^r',extra_dc_swg_h_H_w0p6_Data(1,1),fitted_material_H_swg_h(1,1),'^b')
% % % plot(extra_dc_swg_h_H_w0p6_Data(2,1),sim_neff_H_swg_h(2,1),'^r',extra_dc_swg_h_H_w0p6_Data(2,1),fitted_material_H_swg_h(2,1),'^b')
% % % plot(extra_dc_swg_h_H_w0p6_Data(3,1),sim_neff_H_swg_h(3,1),'^r',extra_dc_swg_h_H_w0p6_Data(3,1),fitted_material_H_swg_h(3,1),'^b')
% % % plot(extra_dc_swg_h_H_w0p6_Data(4,1),sim_neff_H_swg_h(4,1),'^r',extra_dc_swg_h_H_w0p6_Data(4,1),fitted_material_H_swg_h(4,1),'^b')
% % % plot(extra_dc_swg_h_H_w0p6_Data(5,1),sim_neff_H_swg_h(5,1),'^r',extra_dc_swg_h_H_w0p6_Data(5,1),fitted_material_H_swg_h(5,1),'^b')
% % % title('Mode\_MaterialIndex\_vs\_dc\_swg\_h\_w0p6\_sims');
% % % xlabel('dc swg h');
% % % ylabel('Mode Index, Material Index');
% % % dims={'wg width=0.6 micron'};
% % % text(0.5,3, dims);
% % % legend('fitted neff (input:dc\_swg\_h)','effMaterial for fitted neff (input:dc\_swg\_h)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)','sim neff for fitted dc swg h(input is effMaterial))','effmaterial for fitted dc swg h(input is effMaterial))')

% SWG_H  W0P7   
disp('PLOTING Mode and Material vs DC SWG_H W0.7') % IMPORTANT GRAPH
figure('name','SMT_fitted_mode_material_vs_dc_swg_h_w0p7');
plot(dc_swg_h_w0p7,fit_mode_vs_dc_swg_h_w0p7(dc_swg_h_w0p7),'-.r',dc_swg_h_w0p7,fit_material_vs_mode_w0p7(fit_mode_vs_dc_swg_h_w0p7(dc_swg_h_w0p7)),'-.b')
hold on
plot(dc_swg_h_w0p7,neff_swg_h_w0p7,'--c',dc_swg_h_w0p7,fit_material_vs_mode_w0p7(neff_swg_h_w0p7),'--m')
plot(extra_dc_swg_h_H_w0p7_Data(1,1),sim_neff_swg_h_L(1,1),'^r',extra_dc_swg_h_H_w0p7_Data(1,1),fitted_material_swg_h_L_w0p7(1,1),'^b')
plot(extra_dc_swg_h_H_w0p7_Data(2,1),sim_neff_swg_h_L(2,1),'^r',extra_dc_swg_h_H_w0p7_Data(2,1),fitted_material_swg_h_L_w0p7(2,1),'^b')
plot(extra_dc_swg_h_H_w0p7_Data(3,1),sim_neff_swg_h_L(3,1),'^r',extra_dc_swg_h_H_w0p7_Data(3,1),fitted_material_swg_h_L_w0p7(3,1),'^b')
plot(extra_dc_swg_h_H_w0p7_Data(4,1),sim_neff_swg_h_L(4,1),'^r',extra_dc_swg_h_H_w0p7_Data(4,1),fitted_material_swg_h_L_w0p7(4,1),'^b')
plot(extra_dc_swg_h_H_w0p7_Data(5,1),sim_neff_swg_h_L(5,1),'^r',extra_dc_swg_h_H_w0p7_Data(5,1),fitted_material_swg_h_L_w0p7(5,1),'^b')
plot(extra_dc_swg_h_H_w0p7_Data(6,1),sim_neff_swg_h_L(6,1),'^r',extra_dc_swg_h_H_w0p7_Data(6,1),fitted_material_swg_h_L_w0p7(6,1),'^b')
plot(fitted_dc_swg_h_H_vs_mode_w0p7(1,1),sim_neff_H_w0p7(1,1),'^g',fitted_dc_swg_h_H_vs_mode_w0p7(1,1),material_H_w0p7(1,1),'^g')
plot(fitted_dc_swg_h_H_vs_mode_w0p7(2,1),sim_neff_H_w0p7(2,1),'^g',fitted_dc_swg_h_H_vs_mode_w0p7(2,1),material_H_w0p7(2,1),'^g')
plot(fitted_dc_swg_h_H_vs_mode_w0p7(3,1),sim_neff_H_w0p7(3,1),'^g',fitted_dc_swg_h_H_vs_mode_w0p7(3,1),material_H_w0p7(3,1),'^g')
plot(fitted_dc_swg_h_H_vs_mode_w0p7(4,1),sim_neff_H_w0p7(4,1),'^g',fitted_dc_swg_h_H_vs_mode_w0p7(4,1),material_H_w0p7(4,1),'^g')
plot(fitted_dc_swg_h_H_vs_mode_w0p7(5,1),sim_neff_H_w0p7(5,1),'^g',fitted_dc_swg_h_H_vs_mode_w0p7(5,1),material_H_w0p7(5,1),'^g')
plot(fitted_dc_swg_h_H_vs_mode_w0p7(6,1),sim_neff_H_w0p7(6,1),'^g',fitted_dc_swg_h_H_vs_mode_w0p7(6,1),material_H_w0p7(6,1),'^g')


title('Mode\_MaterialIndex\_vs\_dc\_swg\_h\_w0p7\_sims');
xlabel('dc swg h');
ylabel('Mode Index, Material Index');
dims={'wg width=0.7 micron'};
text(0.5,3, dims);
legend('fitted neff (input:dc\_swg\_h)','effMaterial for fitted neff (input:dc\_swg\_h)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)')






% % % % % % % %  REDUNDANT   %%%%%%%%%%%%%%%% TESTS %%%%%%%%
% scenario 1
% SWG_V
% Choose Material Index
% find simulated neff
% find fitted neff
% find fitted dc_swg
% run 3d fdtd with the fitted dc_swg and get the simulated neff ????? need to completed
% % % % disp('PLOTING Mode and Material vs DC SWG_V W0.5')
% % % % figure('name','neff_vs_dc_swg_v_sim_fitted_w0p5');
% % % % plot(dc_swg_v_w0p5,neff_swg_v_w0p5,'r',dc_swg_v_w0p5,fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5),'b');
% % % % hold on
% % % % title('ModeIndex\_dc\_swg\_v\_w0p5\_sims');
% % % % xlabel('dc swg v');
% % % % ylabel('Mode Index (neff)');
% % % % dims={'wg width=0.5 micron'};
% % % % text(0.5,3, dims);
% % % % sim_material_test_value=f_w0p5(:,1);
% % % % sim_neff_test_value=f_w0p5(:,2);
% % % % fitted_neff_test_value=fit_mode_vs_material_w0p5(sim_material_test_value);
% % % % fitted_dc_swg_v_test_value=fit_dc_swg_v_vs_mode_w0p5(sim_neff_test_value);
% % % % plot(fitted_dc_swg_v_test_value,sim_neff_test_value,'o','MarkerEdgeColor','black','MarkerSize',10)
% % % % plot(fitted_dc_swg_v_test_value,fitted_neff_test_value,'p','MarkerEdgeColor','green','MarkerSize',10)
% % % % legend('sim\_neff\_vs\_dc\_swg\_v', 'fitted\_neff\_vs\_dc\_swg\_v','material(input)\_neff(sim)\_dc\_swg\_v(fitted)','material(input)\_neff(fitted)\_dc\_swg\_v(fitted)')
%
%
% % % % % SWG_H
% % % % % Choose Material Index
% % % % % find simulated neff
% % % % % find fitted neff
% % % % % find fitted dc_swg
% % % % % run 3d fdtd with the fitted dc_swg_h and get the simulated neff ????? need to completed
% % % % figure('name','neff_vs_dc_swg_h_sim_fitted_w0p7');
% % % % plot(dc_swg_h_w0p7,neff_swg_h_w0p7,'r',dc_swg_h_w0p7,fit_mode_vs_dc_swg_h_w0p7(dc_swg_h_w0p7),'b');
% % % % hold on
% % % % title('ModeIndex\_dc\_swg\_h\_w0p7\_sims');
% % % % xlabel('dc swg h');
% % % % ylabel('Mode Index (neff)');
% % % % dims={'wg width=0.5 micron'};
% % % % text(0.5,3, dims);
% % % % id1=3;
% % % % sim_material_test_value=f_w0p7(:,1);
% % % % sim_neff_test_value=f_w0p7(:,2);
% % % % fitted_neff_test_value=fit_mode_vs_material_w0p7(sim_material_test_value);
% % % % fitted_dc_swg_h_test_value=fit_dc_swg_h_vs_mode_w0p7(B3,sim_neff_test_value);
% % % % plot(fitted_dc_swg_h_test_value,sim_neff_test_value,'o','MarkerEdgeColor','black','MarkerSize',10)
% % % % plot(fitted_dc_swg_h_test_value,fitted_neff_test_value,'p','MarkerEdgeColor','green','MarkerSize',10)
% % % % legend('sim\_neff\_vs\_dc\_swg\_h', 'fitted\_neff\_vs\_dc\_swg\_h','material(input)\_neff(sim)\_dc\_swg\_h(fitted)','material(input)\_neff(fitted)\_dc\_swg\_h(fitted)')





% scenario 2
% Choose dc_swg_v
% find simulated neff
% find fitted neff
% find fitted Material Index
% run 3d fdtd with the fitted Material Index and get the simulated neff ????? need to completed
figure('name','SMT_neff_vs_material_sim_fitted_w0p6_swg_v'); % IMPORTANT GRAPH
plot(f_w0p6(:,1),f_w0p6(:,2),'-.r',f_w0p6(:,1),fit_mode_vs_material_w0p6(f_w0p6(:,1)),'--b');
hold on
title('ModeIndex\_MaterialIndex\_w0p6\_sims\_swg\_v');
xlabel('Material Index');
ylabel('Mode Index');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);
% id1=3;

plot(material_H(1,1),sim_neff_H(1,1),'p','MarkerEdgeColor','green','MarkerSize',10)
plot(material_H(4,1),sim_neff_H(4,1),'p','MarkerEdgeColor','green','MarkerSize',10)
plot(material_H(5,1),sim_neff_H(5,1),'p','MarkerEdgeColor','green','MarkerSize',10)
plot(fitted_material_swg_v_L_w0p6(1,1),sim_neff_swg_v_L(1,1),'o','MarkerEdgeColor','black','MarkerSize',10)
plot(fitted_material_swg_v_L_w0p6(4,1),sim_neff_swg_v_L(4,1),'o','MarkerEdgeColor','black','MarkerSize',10)
plot(fitted_material_swg_v_L_w0p6(5,1),sim_neff_swg_v_L(5,1),'o','MarkerEdgeColor','black','MarkerSize',10)


% sim_dc_swg_v_test_value=swg_v_data_w0p6(:,1);
% sim_neff_test_value_2=swg_v_data_w0p6(:,2);
% fitted_neff_test_value_2=fit_mode_vs_dc_swg_v_w0p6(sim_dc_swg_v_test_value)
% fitted_material_test_value_2_1=fit_material_vs_mode_w0p6(sim_neff_test_value_2)
% fitted_material_test_value_2_2=fit_material_vs_mode_w0p6(fitted_neff_test_value_2)
% plot(fitted_material_test_value_2_1,sim_neff_test_value_2,'o','MarkerEdgeColor','black','MarkerSize',10)
% plot(fitted_material_test_value_2_2,fitted_neff_test_value_2,'p','MarkerEdgeColor','green','MarkerSize',10)
% error_swg_v_w0p6=100*abs(fitted_material_test_value_2_2-fitted_material_test_value_2_1)./fitted_material_test_value_2_1
legend('sim\_neff\_vs\_material\_index', 'fitted\_neff\_vs\_material\_index')



% % % % figure('name','neff_vs_material_sim_fitted_w0p5_swg_v');
% % % % plot(f_w0p5(:,1),f_w0p5(:,2),'^r',f_w0p5(:,1),fit_mode_vs_material_w0p5(f_w0p5(:,1)),'b');
% % % % hold on
% % % % title('ModeIndex\_MaterialIndex\_w0p5\_sims\_swg\_v');
% % % % xlabel('Material Index');
% % % % ylabel('Mode Index');
% % % % dims={'wg width=0.5 micron'};
% % % % text(0.5,3, dims);
% % % % % id1=3;
% % % % sim_dc_swg_v_test_value=swg_v_data_w0p5(:,1);
% % % % sim_neff_test_value_2=swg_v_data_w0p5(:,2);
% % % % fitted_neff_test_value_2=fit_mode_vs_dc_swg_v_w0p5(sim_dc_swg_v_test_value)
% % % % fitted_material_test_value_2_1=fit_material_vs_mode_w0p5(sim_neff_test_value_2)
% % % % fitted_material_test_value_2_2=fit_material_vs_mode_w0p5(fitted_neff_test_value_2)
% % % % plot(fitted_material_test_value_2_1,sim_neff_test_value_2,'o','MarkerEdgeColor','black','MarkerSize',10)
% % % % plot(fitted_material_test_value_2_2,fitted_neff_test_value_2,'p','MarkerEdgeColor','green','MarkerSize',10)
% % % % legend('sim\_neff\_vs\_material\_index', 'fitted\_neff\_vs\_material\_index','material(fitted)\_neff(sim)\_dc\_swg\_v(input)','material(fitted)\_neff(fitted)\_dc\_swg\_v(input)')
% % % % error_swg_v_w0p5=100*abs(fitted_material_test_value_2_2-fitted_material_test_value_2_1)./fitted_material_test_value_2_1
% % % % 
% % % % figure('name','neff_vs_material_sim_fitted_w0p5_swg_h');
% % % % plot(f_w0p5(:,1),f_w0p5(:,2),'^r',f_w0p5(:,1),fit_mode_vs_material_w0p5(f_w0p5(:,1)),'b');
% % % % hold on
% % % % title('ModeIndex\_MaterialIndex\_w0p5\_sims\_swg\_h');
% % % % xlabel('Material Index');
% % % % ylabel('Mode Index');
% % % % dims={'wg width=0.5 micron'};
% % % % text(0.5,3, dims);
% % % % % id1=3;
% % % % sim_dc_swg_h_test_value=swg_h_data_w0p5([2,10,17,23,25,27],1);
% % % % sim_neff_test_value_3=swg_h_data_w0p5([2,10,17,23,25,27],2);
% % % % fitted_neff_test_value_3=fit_mode_vs_dc_swg_h_w0p5(sim_dc_swg_h_test_value);
% % % % fitted_material_test_value_3_1=fit_material_vs_mode_w0p5(sim_neff_test_value_3);
% % % % fitted_material_test_value_3_2=fit_material_vs_mode_w0p5(fitted_neff_test_value_3);
% % % % error=100*abs(fitted_material_test_value_3_2-fitted_material_test_value_3_1)./fitted_material_test_value_3_1
% % % % 
% % % % plot(fitted_material_test_value_3_1,sim_neff_test_value_3,'o','MarkerEdgeColor','black','MarkerSize',10)
% % % % plot(fitted_material_test_value_3_2,fitted_neff_test_value_3,'p','MarkerEdgeColor','green','MarkerSize',10)
% % % % legend('sim\_neff\_vs\_material\_index', 'fitted\_neff\_vs\_material\_index','material(fitted)\_neff(sim)\_dc\_swg\_h(input)','material(fitted)\_neff(fitted)\_dc\_swg\_h(input)')


figure('name','SMT_neff_vs_material_sim_fitted_w0p7_swg_h');  % IMPORTANT GRAPH
plot(f_w0p7(:,1),f_w0p7(:,2),'-.r',f_w0p7(:,1),fit_mode_vs_material_w0p7(f_w0p7(:,1)),'--b');
hold on
title('ModeIndex\_MaterialIndex\_w0p7\_sims\_swg\_h');
xlabel('Material Index');
ylabel('Mode Index');
dims={'wg width=0.7 micron'};
text(0.5,3, dims);
% id1=3;
plot(material_H_w0p7(3,1),sim_neff_H_swg_h_w0p7(3,1),'p','MarkerEdgeColor','green','MarkerSize',10)
plot(material_H_w0p7(4,1),sim_neff_H_swg_h_w0p7(4,1),'p','MarkerEdgeColor','green','MarkerSize',10)
plot(material_H_w0p7(5,1),sim_neff_H_swg_h_w0p7(5,1),'p','MarkerEdgeColor','green','MarkerSize',10)
plot(fitted_material_swg_h_L_w0p7(3,1),sim_neff_swg_h_L(3,1),'o','MarkerEdgeColor','black','MarkerSize',10)
plot(fitted_material_swg_h_L_w0p7(4,1),sim_neff_swg_h_L(4,1),'o','MarkerEdgeColor','black','MarkerSize',10)
plot(fitted_material_swg_h_L_w0p7(5,1),sim_neff_swg_h_L(5,1),'o','MarkerEdgeColor','black','MarkerSize',10)

% sim_dc_swg_h_test_value=swg_h_data_w0p7(5:end,1);
% sim_neff_test_value_4=swg_h_data_w0p7(5:end,2);
% fitted_neff_test_value_4=fit_mode_vs_dc_swg_h_w0p7(sim_dc_swg_h_test_value);
% fitted_material_test_value_4_1=fit_material_vs_mode_w0p7(sim_neff_test_value_4);
% fitted_material_test_value_4_2=fit_material_vs_mode_w0p7(fitted_neff_test_value_4);
% error_swg_h_w0p7=100*abs(fitted_material_test_value_4_2-fitted_material_test_value_4_1)./fitted_material_test_value_4_1
% 
% plot(fitted_material_test_value_4_1,sim_neff_test_value_4,'o','MarkerEdgeColor','black','MarkerSize',10)
% plot(fitted_material_test_value_4_2,fitted_neff_test_value_4,'p','MarkerEdgeColor','green','MarkerSize',10)
legend('sim\_neff\_vs\_material\_index', 'fitted\_neff\_vs\_material\_index')






