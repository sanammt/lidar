%%%%%%%% Start of Copyright Notice %%%%%%%
%
% copyright © 2018
% Written by Sanam Moslemi-Tabrizi
% All rights reserved
% Can NOT be copied or distrutated without including this copyright notice
%
%%%%%%%% End of Copyright Notice %%%%%%%%%%%
pitch=0.6;
lam = 1.55;
neff = 3.43;
ncl=1.45
factor = lam/neff
target_length=15;
(pitch*neff-lam)/pitch
theta = asind((pitch*neff-lam)/ncl*pitch)
theta_2 = 90-theta
n_periods = ceil(target_length/(pitch))
n_periods_2 = ceil(n_periods/2)
duty_cycle = 0.7;
pitch_2 = 2*pitch
fill_1 = pitch_2 * 0.6 *duty_cycle
etch_1=pitch_2 * 0.6 *(1-duty_cycle)
fill_2 = pitch_2 * 0.4 *duty_cycle
etch_2=pitch_2 * 0.4 *(1-duty_cycle)
L_2 = n_periods_2 * pitch_2 +etch_1
L_2 = n_periods * pitch +etch_1
for i=1:10
    x_min_1(i)=pitch_2*(i-1)+etch_1;
    x_max_1(i)=pitch_2*(i-1)+pitch_2*0.6;
    x_min_2(i)=pitch_2*(i-1)+pitch_2*0.6+etch_2;
    x_max_2(i)=pitch_2*(i);
end

duty_cycle_1 = 0.9;
duty_cycle_2 = 0.8;
duty_cycle_3 = 0.7;
duty_cycle_4 = 0.6;
duty_cycle_5 = 0.5;
fill_1 = pitch * duty_cycle1;
etch_1=pitch *(1-duty_cycle1);
fill_2 = pitch * duty_cycle2;
etch_2=pitch *(1-duty_cycle2);
fill_3 = pitch * duty_cycle3;
etch_3=pitch *(1-duty_cycle3);
fill_4 = pitch * duty_cycle4;
etch_4=pitch *(1-duty_cycle4);
fill_5 = pitch * duty_cycle5;
etch_5=pitch *(1-duty_cycle5);

% non uniform groove size
    x_min_1=etch_1;
    x_max_1=pitch;
    x_min_2=pitch+etch_2;
    x_max_2=pitch*2;
    x_min_3=pitch*2+etch_3;
    x_max_3=pitch*3;
    x_min_4=pitch*3+etch_4;
    x_max_4=pitch*4;
    x_min_5=pitch*4+etch_5;
    x_max_5=pitch*5;
    
    for i=6:n_periods
        x_min(i)=pitch_*(i-1)+etch_5;
        x_max(i)=pitch*(i);
    end
        



