close all
clear all
runDCswg='runDCswg_v_2d.txt';
runDCswgData=load(runDCswg);
figure('name','neff_vs_DCswg_Pitchswg0p1_W0p4_h0p22');
plot(runDCswgData(:,1),runDCswgData(:,2),runDCswgData(:,1),runDCswgData(:,3));
legend('Mode Solution','FDTD 2D');
xlabel('swg DC, pitch=0.1u')
ylabel('Index of Mode, propagation')