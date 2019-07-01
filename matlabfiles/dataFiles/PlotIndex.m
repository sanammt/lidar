close all
clear all

% indexDataFile='IndexVSdc.txt';
% mat2modeDataFile='material2mode.txt';
% %
% indexData=load(indexDataFile);
% figure('name','Index_vs_DC_pitch_100nm_2D_fdtd_res4');
% plot(indexData(:,1),indexData(:,2));
% %
% mat2modeDataData=load(mat2modeDataFile);
% figure('name','material_vs_mode_width_400nm_2D_fdtd_res4');
% plot(mat2modeDataData(:,1),mat2modeDataData(:,2));
% xlabel('Index of Material')
% ylabel('Index of Mode, propagation')
% %%%%%%%
% Mode solutions results 
% file: CheckModes_v.lms
% sweep: neff_Pithswg0p1_W0p4_vs_DCswg_Lswg
% sweep: neff_Pithswg0p1_W0p4_vs_Lswg_DCswg
DC=0:0.05:1;
L=0.5:0.1:1.5;
sizeDC=size(DC);
sizeL=size(L);
neff_vs_DCswg_Lswg='neff_Pitch0p1_W0p4_DC_L_swg.txt';
neff_vs_Lswg_DCswg='neff_Pitch0p1_W0p4_L_DC_swg.txt';
neff_vs_DCswg_LswgData=load(neff_vs_DCswg_Lswg);
neff_vs_Lswg_DCswgData=load(neff_vs_Lswg_DCswg);
%CM = jet(21); 
C = {'k','b','r','g','y','c','m',[.5 .2 .2],[.5 .7 .7],[.8 .2 .6],[0 0.6 .4]};
MK={'+','o','<','^','v','d','x','s','>','*','.'};
figure('name','neff_vs_DCswg_Lswg_Pitchswg0p1_W0p4_h0p22');
plot(L'',neff_vs_DCswg_LswgData(:,1),'LineWidth',1,'color',C{1},'marker','+');
xlabel('Length of grating')
ylabel('Index of Mode, propagation')
hold
for i=2:sizeDC(2)
    %plot(L',neff_vs_DCswg_LswgData(:,i),'color',CM(i,:));
    if (i<=11)
        plot(L',neff_vs_DCswg_LswgData(:,i),'LineWidth',1,'color',C{i},'marker','+');
    elseif (i>11 & i<=22)
        plot(L',neff_vs_DCswg_LswgData(:,i),'LineWidth',1,'color',C{i-11});
    else
        plot(L',neff_vs_DCswg_LswgData(:,i),'LineWidth',1,'color',C{i-11},'marker','o');
    end
end
title('neff\_vs\_DCswg\_Lswg\_Pitchswg0p1\_W0p\_h0p22');
Legend=cell(21,1);
 for iter=1:21
   Legend{iter}=strcat('DC=', num2str(DC(iter)));
 end
 legend(Legend,'location','north');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure('name','neff_vs_Lswg_DCswg_Pitchswg0p1_W0p4_h0p22');
plot(DC',neff_vs_Lswg_DCswgData(:,1),'LineWidth',1,'color',C{1},'marker',MK{1});
xlabel('Duty Cycle (DC)')
ylabel('Index of Mode, propagation')
% legendText=sprintf('DC=%2.2f',DC(1));
% legend(legendText)
hold
for i=2:sizeL(2)
%     plot(L',neff_vs_DCswg_LswgData(:,i),'color',CM(i,:));
    if (i<=5)
        plot(DC',neff_vs_Lswg_DCswgData(:,i),'LineWidth',1,'color',C{i},'marker',MK{i});
    elseif (i>5 & i<=10)
        plot(DC',neff_vs_Lswg_DCswgData(:,i),'LineWidth',1,'color',C{i});
    else
        plot(DC',neff_vs_Lswg_DCswgData(:,i),'LineWidth',1,'color',C{i});
    end
end
title('neff\_vs\_Lswg\_DCswg\_Pitchswg0p1\_W0p4\_h0p22');
Legend=cell(11,1);
 for iter=1:11
   Legend{iter}=strcat('L=', num2str(L(iter)));
 end
 legend(Legend,'location','northwest');

