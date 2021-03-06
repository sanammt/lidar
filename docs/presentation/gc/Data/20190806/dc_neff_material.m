%width_modeSolution
% neff
% material index
close all;
clear all;

materail_neff= cell(50,50);
materail_neff_mod= cell(50,50);

% Read simulation data of 3d fdtd swg_v run.
% data format:  dc_swg - neff - width - pitch_swg
swg_v_data_w0p5=load("fdtd_mode_dc_swg_multi_w0p5_Data.txt");
swg_v_data_w0p6=load("fdtd_mode_dc_swg_multi_w0p6_Data.txt");
%swg_v_data_w0p7=load("fdtd_mode_dc_swg_multi_w0p5_Data.txt");
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




%%%%%%%%%%%%%% CURVE FITTING %%%%%%%%%%%%%

% SWG_V curve fitting
% for w=0.5 ==> id=2 SWG_V
fit_mode_vs_dc_swg_v_w0p5=fit(swg_v_data_w0p5(:,1),swg_v_data_w0p5(:,2),'poly6'); % To get fitted neff for given dc swg_v
fit_dc_swg_v_vs_mode_w0p5=fit(swg_v_data_w0p5(:,2),swg_v_data_w0p5(:,1),'poly6'); % To get fitted dc swg_v for given neff
figure('name','mode_vs_dc_swg_v_w0p5_sim_fitted');
plot(dc_swg_v_w0p5,neff_swg_v_w0p5);
hold on
plot(dc_swg_v_w0p5,fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_v\_w0p5\_sim\_fitted');
xlabel('swg duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
%
% for w=0.6 ==> id=4 SWG_V
fit_mode_vs_dc_swg_v_w0p6=fit(swg_v_data_w0p6(:,1),swg_v_data_w0p6(:,2),'poly6'); % To get fitted neff for given dc swg_v
fit_dc_swg_v_vs_mode_w0p6=fit(swg_v_data_w0p6(:,2),swg_v_data_w0p6(:,1),'poly6'); % To get fitted dc swg_v for given neff
figure('name','mode_vs_dc_swg_v_w0p6_sim_fitted');
plot(dc_swg_v_w0p6,neff_swg_v_w0p6);
hold on
plot(dc_swg_v_w0p6,fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6))
legend('sim','fitted')
title('mode\_vs\_dc\_swg\_v\_w0p6\_sim\_fitted');
xlabel('swg duty cycle');
ylabel('Neff(mode index), sim and fitted');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);



% 3D FDTD curve fitting
% for w=0.5 ==> id=2
f_w0p5=materail_neff_mod{2};
fit_mode_vs_material_w0p5=fit(f_w0p5(:,1),f_w0p5(:,2),'linearinterp'); % To get fitted neff for given effMaterial Index
fit_material_vs_mode_w0p5=fit(f_w0p5(:,2),f_w0p5(:,1),'linearinterp'); % To get fitted effMaterial Index for given neff 
%
% for w=0.6 ==> id=4
f_w0p6=materail_neff_mod{4};
fit_mode_vs_material_w0p6=fit(f_w0p6(:,1),f_w0p6(:,2),'linearinterp'); % To get fitted neff for given effMaterial Index
fit_material_vs_mode_w0p6=fit(f_w0p6(:,2),f_w0p6(:,1),'linearinterp'); % To get fitted effMaterial Index for given neff 

% For a given width, for given values of dc_swg_v, we have neff obtained from 3d fdtd swg_v simulation
%********** dc_swg ==> sims: neff ==> fitted 3d FDTD: fitted effMaterialIndex
% which means for given dc_swg we can calculated the respective effMaterialIndex (for a fixed defined width)
% neff_swg_v_w0p5 is the simulated data for dc_swg sweep
% Hybrid mode: mix of simulation and curve fitting
material_vs_dc_swg_w0p5=fit_material_vs_mode_w0p5(neff_swg_v_w0p5);
material_vs_dc_swg_w0p6=fit_material_vs_mode_w0p6(neff_swg_v_w0p6);
%
figure('name','mode_material_vs_dc_swg_v_w0p5_hybrid_sims_fitted');
plot(dc_swg_v_w0p5,neff_swg_v_w0p5,'r',dc_swg_v_w0p5,material_vs_dc_swg_w0p5,'b');
legend('Mode Index','Material Index');
title('ModeIndex\_MaterialIndex\_dc\_swg\_v\_w0p5\_hybrid\_sims\_fitted');
xlabel('swg duty cycle');
ylabel('Mode Index , Material Index');
dims={'wg width=0.5 micron'};
text(0.5,3, dims);
% 
% For a given width, for given values of dc_swg_v, we have neff obtained from 3d fdtd swg_v simulation
%********** dc_swg ==> sims: neff ==> fitted 3d FDTD: fitted effMaterialIndex
% which means for given dc_swg we can calculated the respective effMaterialIndex (for a fixed defined width)
% neff_swg_v_w0p6 is the simulated data for dc_swg sweep
% Hybrid mode: mix of simulation and curve fitting
material_vs_dc_swg_w0p6=fit_material_vs_mode_w0p6(neff_swg_v_w0p6);
%
figure('name','mode_material_vs_dc_swg_v_w0p6_hybrid_sims_fitted');
plot(dc_swg_v_w0p6,neff_swg_v_w0p6,'r',dc_swg_v_w0p6,material_vs_dc_swg_w0p6,'b');
legend('Mode Index','Material Index');
title('ModeIndex\_MaterialIndex\_dc\_swg\_v\_w0p6\_hybrid\_sims\_fitted');
xlabel('swg duty cycle');
ylabel('Mode Index , Material Index');
dims={'wg width=0.6 micron'};
text(0.5,3, dims);
% 
% %
% 
% % for w=0.7 ==> id=6
% f_w0p7=materail_neff_mod{6};
% fit_mode_vs_material_w0p7=fit(f_w0p7(:,1),f_w0p7(:,2),'linearinterp');
% fit_material_vs_mode_w0p7=fit(f_w0p7(:,2),f_w0p7(:,1),'linearinterp');
% 
% material_vs_dc_swg_w0p7=fit_material_vs_mode_w0p7(neff_swg_v_w0p7);
% figure('name','mode_material_dc_swg_v_w0p7');
% plot(dc_swg_v_w0p7,neff_swg_v_w0p7,'r',dc_swg_v_w0p7,material_vs_dc_swg_w0p7,'b');
% legend('Mode Index','Material Index');
% title('ModeIndex\_MaterialIndex\_dc\_swg\_v\_w0p7');
% xlabel('swg duty cycle');
% ylabel('Mode Index , Material Index');
% dims={'wg width=0.7 micron'};
% text(0.5,3, dims);


% FULL CURVE FITTING
% w0p5
% Input is dc_swg_v => Calculate the fitted neff => Calculate the fitted effMaterialIndex
neff_fitted_w0p5_1=fit_mode_vs_dc_swg_v_w0p5(dc_swg_v_w0p5);
material_fitted_w0p5_1=fit_material_vs_mode_w0p5(neff_fitted_w0p5_1);

figure('name','fitted_mode_material_vs_dc_swg_v_w0p5');
%'fitted neff','effMaterial for fitted neff'
plot(dc_swg_v_w0p5,neff_fitted_w0p5_1,'r',dc_swg_v_w0p5,material_fitted_w0p5_1,'b')
hold on
%'sim neff','fitted effmaterial for sim neff'
plot(dc_swg_v_w0p5,neff_swg_v_w0p5,'*r',dc_swg_v_w0p5,material_vs_dc_swg_w0p5,'*b')
%legend('fitted neff (input:dc\_swg\_v)','effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)')


% Input is effMaterial => Calculate the fitted neff => Calculate the fitted dc_swg_v
neff_fitted_w0p5_2=fit_mode_vs_material_w0p5(material_vs_dc_swg_w0p5);
dc_swg_v_w0p5_fitted_2=fit_dc_swg_v_vs_mode_w0p5(neff_fitted_w0p5_2);
plot(dc_swg_v_w0p5_fitted_2,neff_fitted_w0p5_2,'^r',dc_swg_v_w0p5_fitted_2,material_vs_dc_swg_w0p5,'^b')
legend('fitted neff (input:dc\_swg\_v)','effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)','fitted neff (input:given effMaterial Index)','given effMaterial Index')


% error_dc_swg=100*(abs(dc_swg_v_w0p5_fitted_2-dc_swg_v_w0p5)./dc_swg_v_w0p5)

% figure('name','fitted_mode_material_vs_dc_swg_v_w0p5');
% plot(dc_swg_v_w0p5_fitted_2,neff_fitted_2,'r',dc_swg_v_w0p5_fitted_2,material_vs_dc_swg_w0p5,'b')
% hold on
% plot(dc_swg_v_w0p5_fitted_2,neff_swg_v_w0p5,'*r',dc_swg_v_w0p5,material_vs_dc_swg_w0p5,'*b')
% legend('fitted neff','material for fitted neff','sim neff','fitted material for sim neff')


% w0p6
% Input is dc_swg_v => Calculate the fitted neff => Calculate the fitted effMaterialIndex
neff_fitted_w0p6_1=fit_mode_vs_dc_swg_v_w0p6(dc_swg_v_w0p6);
material_fitted_w0p6_1=fit_material_vs_mode_w0p6(neff_fitted_w0p6_1);

figure('name','fitted_mode_material_vs_dc_swg_v_w0p6');
%'fitted neff','effMaterial for fitted neff'
plot(dc_swg_v_w0p6,neff_fitted_w0p6_1,'r',dc_swg_v_w0p6,material_fitted_w0p6_1,'b')
hold on
%'sim neff','fitted effmaterial for sim neff'
plot(dc_swg_v_w0p6,neff_swg_v_w0p6,'*r',dc_swg_v_w0p6,material_vs_dc_swg_w0p6,'*b')
%legend('fitted neff (input:dc\_swg\_v)','effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)')


% Input is effMaterial => Calculate the fitted neff => Calculate the fitted dc_swg_v
neff_fitted_w0p6_2=fit_mode_vs_material_w0p6(material_vs_dc_swg_w0p6);
dc_swg_v_w0p6_fitted_2=fit_dc_swg_v_vs_mode_w0p6(neff_fitted_w0p6_2);
plot(dc_swg_v_w0p6_fitted_2,neff_fitted_w0p6_2,'^r',dc_swg_v_w0p6_fitted_2,material_vs_dc_swg_w0p6,'^b')
% legend('fitted neff (input:dc\_swg\_v)','effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)','fitted neff (input:given effMaterial Index)','given effMaterial Index')


% error_dc_swg_w0p5=100*(abs(dc_swg_v_w0p5_fitted_2-dc_swg_v_w0p5)./dc_swg_v_w0p5);
% error_dc_swg_w0p6=100*(abs(dc_swg_v_w0p6_fitted_2-dc_swg_v_w0p6)./dc_swg_v_w0p6);

%%%%% EXTRA DATA %%%%%
extra_mat_neff=[2.5342 , 1.58445 ; 2.6224,1.65227 ; 2.6866,1.70555 ; 2.7738,1.7851 ; 2.8702,1.8808 ]
extra_fitted_dc_swg_vs_neff=fit_dc_swg_v_vs_mode_w0p6(extra_mat_neff(:,2))
extra_material_vs_neff=fit_material_vs_mode_w0p6(extra_mat_neff(:,2))
extra_fitted_neff_vs_material=fit_mode_vs_material_w0p6(extra_material_vs_neff)
%%%%%%%%%%%%%%%%%%%%%%
plot(extra_fitted_dc_swg_vs_neff(1),extra_fitted_neff_vs_material(1,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)
plot(extra_fitted_dc_swg_vs_neff(1),extra_material_vs_neff(1,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)

plot(extra_fitted_dc_swg_vs_neff(1),extra_mat_neff(1,1),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)
plot(extra_fitted_dc_swg_vs_neff(1),extra_mat_neff(1,2),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)

%
plot(extra_fitted_dc_swg_vs_neff(2),extra_fitted_neff_vs_material(2,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)
plot(extra_fitted_dc_swg_vs_neff(2),extra_material_vs_neff(2,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)

plot(extra_fitted_dc_swg_vs_neff(2),extra_mat_neff(2,1),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)
plot(extra_fitted_dc_swg_vs_neff(2),extra_mat_neff(2,2),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)

%
plot(extra_fitted_dc_swg_vs_neff(3),extra_fitted_neff_vs_material(3,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)
plot(extra_fitted_dc_swg_vs_neff(3),extra_material_vs_neff(3,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)

plot(extra_fitted_dc_swg_vs_neff(3),extra_mat_neff(3,1),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)
plot(extra_fitted_dc_swg_vs_neff(3),extra_mat_neff(3,2),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)

%
plot(extra_fitted_dc_swg_vs_neff(4),extra_fitted_neff_vs_material(4,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)
plot(extra_fitted_dc_swg_vs_neff(4),extra_material_vs_neff(4,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)

plot(extra_fitted_dc_swg_vs_neff(4),extra_mat_neff(4,1),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)
plot(extra_fitted_dc_swg_vs_neff(4),extra_mat_neff(4,2),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)

%
plot(extra_fitted_dc_swg_vs_neff(5),extra_fitted_neff_vs_material(5,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)
plot(extra_fitted_dc_swg_vs_neff(5),extra_material_vs_neff(5,1),'-o','MarkerEdgeColor','black','MarkerSize',5,'LineWidth',2)

plot(extra_fitted_dc_swg_vs_neff(5),extra_mat_neff(5,1),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)
plot(extra_fitted_dc_swg_vs_neff(5),extra_mat_neff(5,2),'-x','MarkerEdgeColor','green','MarkerSize',5,'LineWidth',1)

%
legend('fitted neff (input:dc\_swg\_v)','effMaterial for fitted neff (input:dc\_swg\_v)','sim neff (dc sweep sims)','effmaterial for sim neff (dc sweep sims)','fitted neff (input:given effMaterial Index)','given effMaterial Index')
%
%
%
figure('name','extra');
plot(material_vs_dc_swg_w0p6,neff_fitted_w0p6_2)
hold on
plot(extra_material_vs_neff(1),extra_fitted_neff_vs_material(1,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)
plot(extra_mat_neff(1,1),extra_mat_neff(1,2),'-p','MarkerEdgeColor','red','MarkerSize',10,'LineWidth',2)
%
plot(extra_material_vs_neff(2),extra_fitted_neff_vs_material(2,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)
plot(extra_mat_neff(2,1),extra_mat_neff(2,2),'-p','MarkerEdgeColor','red','MarkerSize',10,'LineWidth',2)
%
plot(extra_material_vs_neff(3),extra_fitted_neff_vs_material(3,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)
plot(extra_mat_neff(3,1),extra_mat_neff(3,2),'-p','MarkerEdgeColor','red','MarkerSize',10,'LineWidth',2)
%
plot(extra_material_vs_neff(4),extra_fitted_neff_vs_material(4,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)
plot(extra_mat_neff(4,1),extra_mat_neff(4,2),'-p','MarkerEdgeColor','red','MarkerSize',10,'LineWidth',2)
%
plot(extra_material_vs_neff(5),extra_fitted_neff_vs_material(5,1),'-x','MarkerEdgeColor','black','MarkerSize',10,'LineWidth',2)
plot(extra_mat_neff(5,1),extra_mat_neff(5,2),'-p','MarkerEdgeColor','red','MarkerSize',10,'LineWidth',2)
%


