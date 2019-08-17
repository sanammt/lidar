close all;
clear all;
C = {'k','b','r','g','y','c','m',[.5 .2 .2],[.5 .7 .7],[.8 .2 .6],[0 0.6 .4]};
MK={'+','o','<','^','v','d','x','s','>','*','.'};

w_neff_modeSolution=zeros(6,34);
materail_neff= cell(50,50);

fdtddata=load("fdtd_effMaterial_Index_effMaterial_multi_Data.txt");
fdtddata1=load("fdtd_effMaterial_Index_effMaterial1_multi_Data.txt");
fdtddata2=load("fdtd_effMaterial_Index_effMaterial2_multi_Data.txt");
modedata=load("effMaterial_Index_widerSio2_effMaterial_multi_Data.txt");
modedata1=load("effMaterial_Index_widerSio2_effMaterial1_multi_Data.txt");
modedata2=load("effMaterial_Index_widerSio2_effMaterial2_multi_Data.txt");
modedata3=load("effMaterial_Index_widerSio2_effMaterial3_multi_Data.txt");
modedata4=load("effMaterial_Index_widerSio2_effMaterial4_multi_Data.txt");
modedata5=load("effMaterial_Index_widerSio2_effMaterial5_multi_Data.txt");
%
sio2Index=1.44*ones(34,1);
materialIndex=modedata(:,1);
materialIndex1=modedata1(:,1);
materialIndex2=modedata2(:,1);
materialIndex3=modedata3(:,1);
materialIndex4=modedata4(:,1);
materialIndex5=modedata5(:,1);
%
neff_modeMonitor_fdtd=fdtddata(:,2);
neff_modeMonitor_fdtd1=fdtddata1(:,2);
neff_modeMonitor_fdtd2=fdtddata2(:,2);
%
neff_3d_fdtd=fdtddata(:,3);
neff_3d_fdtd1=fdtddata1(:,3);
neff_3d_fdtd2=fdtddata2(:,3);
%
neff_modeSolution=modedata(:,2);
neff_modeSolution1=modedata1(:,2);
neff_modeSolution2=modedata2(:,2);
neff_modeSolution3=modedata3(:,2);
neff_modeSolution4=modedata4(:,2);
neff_modeSolution5=modedata5(:,2);
%
Nrow0 = min(length(neff_modeMonitor_fdtd),length(neff_modeMonitor_fdtd1));
Nrow1=min(Nrow0,length(neff_modeMonitor_fdtd2))
Nrow=min(Nrow1,length(neff_modeSolution))
width_modeSolution=modedata(:,3);
width_modeSolution1=modedata1(:,3);
width_modeSolution2=modedata2(:,3);
width_modeSolution3=modedata3(:,3);
width_modeSolution4=modedata4(:,3);
width_modeSolution5=modedata5(:,3);

figure('name','ModeIndex_vs_Waveguide_Width_Mode_Solution');
set(gcf, 'Position',  [100, 100, 1000, 1000])
plot(width_modeSolution,neff_modeSolution,'r',width_modeSolution1,neff_modeSolution1,'b',width_modeSolution2,neff_modeSolution2,'g',width_modeSolution3,neff_modeSolution3,'k',width_modeSolution4,neff_modeSolution4,'m',width_modeSolution5,neff_modeSolution5,'c');
% hold on
% plot(width_modeSolution,neff_modeSolution,'.r',width_modeSolution1,neff_modeSolution1,'.b',width_modeSolution2,neff_modeSolution2,'.g',width_modeSolution3,neff_modeSolution3,'.k',width_modeSolution4,neff_modeSolution4,'.m',width_modeSolution5,neff_modeSolution5,'.c');
legend('(Silicon), Mode','(3.25), Mode','(3.0), Mode','(2.75), Mode','(2.5), Mode','(2.25), Mode','location','northwest');
%legend('(Silicon), Mode','(3.25), Mode','(3.0), Mode','(2.75), Mode','(2.5), Mode','(2.25), Mode','(Silicon), 3D FDFD','(3.25), 3D FDFD','(3.0), 3D FDFD','(2.75), 3D FDFD','(2.5), 3D FDFD','(2.25), 3D FDFD','location','northwest');
title('ModeIndex\_vs\_Waveguide\_Width\_Mode\_Solution');
xlabel('Waveguide Width (micron)');
ylabel('Neff, Propagation Mode Index');
% dims={'wg width=0.5 micron','wg depth=0.22 micron', '' ,'sio2 width=1 micron','sio2 depth=1.5 micron'};
% text(materialIndex(end-1),neff_3d_fdtd(end-3), dims);

figure('name','ModeIndex_vs_Waveguide_Width_Different_simulation_methods');
set(gcf, 'Position',  [100, 100, 1000, 1000])
plot(width_modeSolution(2:Nrow,1),neff_modeSolution(2:Nrow,1),'r',width_modeSolution1(2:Nrow,1),neff_modeSolution1(2:Nrow,1),'b',width_modeSolution2(2:Nrow,1),neff_modeSolution2(2:Nrow,1),'g',width_modeSolution3(2:Nrow,1),neff_modeSolution3(2:Nrow,1),'k',width_modeSolution4(2:Nrow,1),neff_modeSolution4(2:Nrow,1),'m',width_modeSolution5(2:Nrow,1),neff_modeSolution5(2:Nrow,1),'c');
hold on
plot(width_modeSolution(2:Nrow,1),neff_3d_fdtd(2:Nrow,1),'*r',width_modeSolution1(2:Nrow,1),neff_3d_fdtd1(2:Nrow,1),'*b',width_modeSolution1(2:Nrow,1),neff_3d_fdtd2(2:Nrow,1),'*g');
plot(width_modeSolution(2:Nrow,1),neff_modeMonitor_fdtd(2:Nrow,1),'^r',width_modeSolution1(2:Nrow,1),neff_modeMonitor_fdtd1(2:Nrow,1),'^b',width_modeSolution1(2:Nrow,1),neff_modeMonitor_fdtd2(2:Nrow,1),'^g');
legend('(Silicon), Mode Solution','(3.25), Mode Solution','(3.0), Mode Solution','(2.75), Mode Solution','(2.5), Mode Solution','(2.25), Mode Solution','(Silicon), 3D FDFD','(3.25), 3D FDFD','(Silicon), FDFD Mode Monitor','(3.25), 3D FDFD Mode Monitor','(3.0), 3D FDFD Mode Monitor','location','northwest');
%legend('(Silicon), Mode','(3.25), Mode','(3.0), Mode','(2.75), Mode','(2.5), Mode','(2.25), Mode','(Silicon), 3D FDFD','(3.25), 3D FDFD','(3.0), 3D FDFD','(2.75), 3D FDFD','(2.5), 3D FDFD','(2.25), 3D FDFD','location','northwest');
title('ModeIndex\_vs\_Waveguide\_Width\_Different\_simulation\_methods');
xlabel('Waveguide Width (micron)');
ylabel('Neff, Propagation Mode Index');


fdfd_3d_vs_modeExpansion=100*abs(neff_modeMonitor_fdtd(2:Nrow,1)-neff_3d_fdtd(2:Nrow,1))./neff_3d_fdtd(2:Nrow,1);
fdfd_3d_vs_modeExpansion1=100*abs(neff_modeMonitor_fdtd1(2:Nrow,1)-neff_3d_fdtd1(2:Nrow,1))./neff_3d_fdtd1(2:Nrow,1);
fdfd_3d_vs_modeExpansion2=100*abs(neff_modeMonitor_fdtd2(2:Nrow,1)-neff_3d_fdtd2(2:Nrow,1))./neff_3d_fdtd2(2:Nrow,1);

for i=1:Nrow
    w_neff_modeSolution(1,i)=neff_modeSolution(i,1);
    w_neff_modeSolution(2,i)=neff_modeSolution1(i,1);
    w_neff_modeSolution(3,i)=neff_modeSolution2(i,1);
    w_neff_modeSolution(4,i)=neff_modeSolution3(i,1);
    w_neff_modeSolution(5,i)=neff_modeSolution4(i,1);
    w_neff_modeSolution(6,i)=neff_modeSolution5(i,1);
end
mat_neff_modeSolution=[width_modeSolution,neff_modeSolution,neff_modeSolution1,neff_modeSolution2,neff_modeSolution3,neff_modeSolution4,neff_modeSolution5];
material_index_range=[materialIndex(1,1), materialIndex1(1,1), materialIndex2(1,1), materialIndex3(1,1), materialIndex4(1,1), materialIndex5(1,1)];

figure('name','ModeIndex_vs_MaterialIndex_Mode_Solution_All_Width');
set(gcf, 'Position',  [100, 100, 1000, 1000])
hold on
for i=1:34 
    plot(material_index_range,mat_neff_modeSolution(i,2:7));
    materail_neff{i}=[material_index_range',mat_neff_modeSolution(i,2:7)'];
end
title('ModeIndex\_vs\_MaterialIndex\_Mode\_Solution\_All\_Width');
xlabel('Material Index');
ylabel('Neff, Propagation Mode Index');
dim1={'wg width=0.4 micron'};
text(material_index_range(1),mat_neff_modeSolution(1,2), dim1);
dim2={'wg width=0.45 micron'};
text(material_index_range(1),mat_neff_modeSolution(2,2), dim2);
dim3={'wg width=0.5 micron'};
text(material_index_range(1),mat_neff_modeSolution(3,2), dim3);
dim4={'wg width=0.55 micron'};
text(material_index_range(1),mat_neff_modeSolution(4,2), dim4);
dim5={'wg width=0.6 micron'};
text(material_index_range(1),mat_neff_modeSolution(5,2), dim5);
dim6={'wg width=0.65 micron'};
text(material_index_range(1),mat_neff_modeSolution(6,2), dim6);
dim7={'wg width=0.7 micron'};
text(material_index_range(1),mat_neff_modeSolution(7,2), dim7);
dim34={'wg width=2 micron'};
text(material_index_range(1),mat_neff_modeSolution(34,2), dim34);




% fitted=0.9329*material_index_range-0.4548;
% 
% figure;
% plot(material_index_range,mat_neff_modeSolution(34,2:7),material_index_range,mat_neff_modeSolution(30,2:7),material_index_range,mat_neff_modeSolution(25,2:7));
% hold on
% plot(material_index_range,fitted,'*r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plot(materialIndex,neff_modeSolution,'r',materialIndex,neff_modeMonitor_fdtd,'b',materialIndex,neff_3d_fdtd','*g');
% legend('Mode Solution','Mode Monitor in FDTD','3D FDTD','location','northwest');
% title('ModeIndex\_vs\_MaterialIndex\_fixedWidth');
% xlabel('Material Index');
% ylabel('Neff, Propagation Mode Index');
% dims={'wg width=0.5 micron','wg depth=0.22 micron', '' ,'sio2 width=1 micron','sio2 depth=1.5 micron'};
% text(materialIndex(end-1),neff_3d_fdtd(end-3), dims);
