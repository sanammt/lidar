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
ilam=0
for lam=[1.35 1.55 1.65] % micron
    ilam=ilam+1;
%lam =1.55 % micron
d = 0.4 % pitch in micron
dc = 0.6 %duty cycle
tg=0.2 % micron
n_cl = 1.44
eps_cl=n_cl^2
n_etch = 1.44
eps_g=n_etch^2
eps_s = 1.44^2 % substrate (must be isolator) permittivity
k0=2*pi/lam
%%%%%%%%%%%%%%%%
xylimit=4;
syms x y
n_air = 1;
n_s=sqrt(eps_s);
f1=y-x-n_cl; %y-x<n_air => f1<0
f2=x-y-n_cl; %x-y>n_air =f2>0
f3=2*x-y-n_s; % 2*x-y>n_s =>f3>0
f4=y-x; %y~=x => f4~=0
f5=y-0.0001*x-n_s; %y>ns => f5>0
ezplot(y-x==n_cl,[1,xylimit,1,xylimit]) % f1<0
hold on
ezplot(x-y==n_cl,[1,xylimit,1,xylimit]) % f2>0
ezplot(2*x-y==n_s,[1,xylimit,1,xylimit]) % f3>0
ezplot(x-y==0,[1,xylimit,1,xylimit])     % f4~=0
ezplot(y-0.0001*x==n_s,[1,xylimit,1,xylimit]) % f5>0
ezplot(0.0001*y+x==lam/d,[1,xylimit,1,xylimit]) % x=lam/d
title([])
xlabel('\lambda/d \rightarrow','FontSize',16);
ylabel('N \uparrow','FontSize',16);
set(get(gca,'ylabel'),'rotation',0);
% set(get(gca, 'ylabel'), 'Position', get(gca, 'ylabel') + [-0.5 0 0]);
ylabh = get(gca,'ylabel');
set(ylabh,'Position',get(ylabh,'Position') + [-0.1 0.1 0]);

%
clear x
x=linspace(1, xylimit);
y11=x+n_cl; %y-x<n_air => f1<0
y21=x-n_cl; %x-y>n_air =f2>0
y12=2*x-n_s; % 2*x-y>n_s =>f3>0
y3=x.*0+n_s; %y>ns => f5>0
y13=min(y11,y12);
y1=max(y13,y3);
y2=max(y21,y3);
y41=x; %y~=x => f4~=0
y4=max(y41,y3);
%f4=y-x; %y~=x => f4~=0


% plot(x, y1, 'r', 'LineWidth', 2);
 ylim([1 xylimit])
% plot(x, y4, 'r', 'LineWidth', 2);
x2 = [x, fliplr(x)];
inBetween1 = [y1, fliplr(y4)];
fill(x2, inBetween1, 'y');


% plot(x, y2, 'b', 'LineWidth', 2);
inBetween2 = [y2, fliplr(y4)];
fill(x2, inBetween2, 'c');
%txt = '$$\bullet \leftarrow 0.25t e^{-0.005t} at t = 300 \int_{0}^{2} x^2\sin(x) dx $$';
%txt2 ='$$ \int_{0}^{2} x^2\sin(x) dx $$'; 
%txt='hello'
txt1 = '$$\leftarrow \frac{\lambda}{d} - N= n_{cl} $$';
text(3.5,2.5,txt1,'Interpreter','latex')

txt2 = '$$\leftarrow N-\frac{\lambda}{d}= n_{cl} $$';
text(2.75,3.74,txt2,'Interpreter','latex')

txt3 = '$$\leftarrow N = \frac{\lambda}{d} $$';
text(3,3,txt3,'Interpreter','latex')

txt4 = '$$\downarrow N = \sqrt\epsilon_s $$';
text(2,1.5,txt4,'Interpreter','latex')

txt5 = '$$\leftarrow \frac{\lambda}{d} - N = \sqrt\epsilon_s $$';
text(2.1,2.75,txt5,'Interpreter','latex')

txt6 = 'Forward';
text(2.4,3,txt6,'FontSize',14,'Color','red','FontWeight','bold')
txt7 = 'Backward';
text(2.4,2.05,txt7,'FontSize',14,'Color','red','FontWeight','bold')

txt7 = '$$\leftarrow \frac{\lambda}{d}$$';
text(lam/d,sqrt(eps_s)+0.1,txt7,'Interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%
% abs(Neff-lam/d)<1 (1=n_cl?)
Neff_max_1=lam/d + n_cl; % Neff < lam/d + n_cl
Neff_min_temp(ilam)=lam/d - n_cl % Neff > lam/d - n_cl
Neff_max_2 = (2*lam)/d - sqrt(eps_s);
Neff_max_temp(ilam) = max(Neff_max_1,Neff_max_2)
end
Neff_min=max(Neff_min_temp)
Neff_max=min(Neff_max_temp)
lam=1.55 % micron

%%% Get Neff
NeffQ='What is the value for Neff?\n';
Neff=input(NeffQ);

eps_r=Neff^2
cond1 = (Neff ~= lam/d)
beta0=Neff*k0;

% Neff > sqrt(eps_s)
cond2 = (Neff > sqrt(eps_s))
leakage_factor =((eps_r - eps_cl)*sin(pi*dc))^2;
cond3 = (Neff > n_cl) % close to n_cl



condf1= (abs(Neff-lam/d)<n_cl); %y-x<n_air => f1<0
condf3=(2*lam/d-Neff>n_s); % 2*x-y>n_s =>f3>0


NeffOkey1=cond1 & cond2 & cond3 & condf1 & condf3;
while (NeffOkey1 == 0)
   disp(Neff_min)
   disp(Neff_max)
   NeffQ='Repeat: What is the value for Neff?\n';
   Neff=input(NeffQ);
   eps_r=Neff^2;
   cond1 = (Neff ~= lam/d);
   % Neff > sqrt(eps_s) 
   cond2 = (Neff > sqrt(eps_s));
   leakage_factor =((eps_r - eps_cl)*sin(pi*dc))^2;
   cond3 = (Neff > n_cl); % close to n_cl
   condf1= (abs(Neff-lam/d)<n_cl); %y-x<n_air => f1<0
   condf3=(2*lam/d-Neff>n_s); % 2*x-y>n_s =>f3>0
   NeffOkey1=cond1 & cond2 & cond3 & condf1 & condf3;
end
 
dcArray=0.55:0.1:0.7;%[0.55 0.6 0.65 0.7];
NeffArray=1.8:0.1:2.8;%[1.8 2.2 2.6 2.8];
dArray=0.4:0.2:0.8;%[0.45 0.5 0.55 0.6 0.7];
id_tg=0;
for i=1:length(dcArray)
    for j=1:length(NeffArray)
        for k=1:length(dArray)
            id_tg=id_tg+1;
            n_g=dcArray(i)*NeffArray(j)+(1-dcArray(i))*n_etch;
            eps_g=n_g^2;
            dc_tg=dcArray(i);
            Neff_tg=NeffArray(j);
            d_tg=dArray(k);
            tg_temp(id_tg)=(0.25*lam)/sqrt(eps_g-(NeffArray(j)-lam/dArray(k))^2);
%             sprintf('dc=%2.2f Neff=%2.2f d=%2.2f => tg=%2.2f',dcArray(i),NeffArray(j),dArray(k),tg)
        end
    end
end
tgMax=min(tg_temp)
if (tg > tgMax)
    tgQ='What is the value for tg?\n';
    tg=input(tgQ);
end
% dc1=0.6:0.05:0.75;
% leakage_factor = ((eps_g-eps_cl)*sin(dc1*pi)).^2
% figure('name','leakage_factor');
% plot(dc1,leakage_factor)
% 
% 
% lam1=1.35:0.05:1.65;
% Wc_fact=sqrt(1-(Neff-(lam1/d).^2)/eps_cl);
% figure('name','wc_factor');
% plot(lam1,Wc_fact)

lam1=1.35:0.05:1.65;
theta1 = asind((Neff-(lam1/d))/n_cl)
fov1=max(theta1)-min(theta1)

theta2=n_cl*theta1;
fov2=max(theta2)-min(theta2)
sprintf('Neff (effective film)=%4.4f',Neff)
sprintf('tg (etch depth)=%4.4f',tg)
sprintf('d (pitch)=%4.4f',d)
sprintf('dc=%4.4f',dc)
sprintf('n_s (isolator below the thin film )=%4.4f',n_s)
sprintf('n_cl (upper cladding and the groove)=%4.4f',n_cl)




% n_g=dc*Neff+(1-dc)*Neff
% eps_g=n_g^2
% tg=(0.25*lam)/sqrt(eps_g-(Neff-lam/d)^2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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