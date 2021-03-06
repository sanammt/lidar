%%%%%%%% Start of Copyright Notice %%%%%%%
%
% copyright � 2018
% Written by Sanam Moslemi-Tabrizi
% All rights reserved
% Can NOT be copied or distrutated without including this copyright notice
%
%%%%%%%% End of Copyright Notice %%%%%%%%%%%


clear all
close all
xylimit=4;
% lam =1.55 % micron
% d = 0.500 % pitch in micron
% dc = 0.6 %duty cycle
% n_cl = 1.44
% eps_cl=sqrt(n_cl)
% eps_s = 1.44^2 
% eps_r=2.6^2
% 
% 
% cond1 = (Neff ~= lam/d)
% 
% % Neff > sqrt(eps_s)
% cond2 = (Neff > sqrt(eps_s))
% leakage_factor =((eps_r - eps_cl)*sin(pi*dc))^2;
% cond3 = (Neff > n_cl) % close to n_cl
% x=1:0.5:5;
% y=1:0.5:5;

%[x,y]=meshgrid(1:0.5:5)
syms x y
% ezplot(x^2==y^4)
n_air = 1;
n_s=1.44
f1=y-x-n_air; %y-x<n_air => f1<0
f2=x-y-n_air; %x-y>n_air =f2>0
f3=2*x-y-n_s; % 2*x-y>n_s =>f3>0
f4=y-x; %y~=x => f4~=0
f5=y-0.0001*x-n_s; %y>ns => f5>0
%
% patch(f1,f2,'r');
ezplot(y-x==n_air,[1,xylimit,1,xylimit]) % f1<0
hold on
ezplot(x-y==n_air,[1,xylimit,1,xylimit]) % f2>0
ezplot(2*x-y==n_s,[1,xylimit,1,xylimit]) % f3>0
ezplot(x-y==0,[1,xylimit,1,xylimit])     % f4~=0
ezplot(y-0.0001*x==n_s,[1,xylimit,1,xylimit]) % f5>0
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

n_air = 1;
n_s=1.44

y11=x+n_air; %y-x<n_air => f1<0
y21=x-n_air; %x-y>n_air =f2>0
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
txt1 = '$$\leftarrow \frac{\lambda}{d} - N= 1 $$';
text(3,2,txt1,'Interpreter','latex')

txt2 = '$$\leftarrow N-\frac{\lambda}{d}= 1 $$';
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



%%%% fill works well

% clear all
% close all
% x=linspace(1, 5);
% 
% n_air = 1;
% n_s=1.44
% 
% y11=x+n_air; %y-x<n_air => f1<0
% y21=x-n_air; %x-y>n_air =f2>0
% y12=2*x-n_s; % 2*x-y>n_s =>f3>0
% y3=x.*0+n_s; %y>ns => f5>0
% y13=min(y11,y12);
% y1=max(y13,y3);
% y2=max(y21,y3);
% %f4=y-x; %y~=x => f4~=0
% 
% 
% plot(x, y1, 'r', 'LineWidth', 2);
% hold on;
% plot(x, y2, 'b', 'LineWidth', 2);
% x2 = [x, fliplr(x)];
% inBetween = [y1, fliplr(y2)];
% fill(x2, inBetween, 'g');
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% patch([x(1:x_y1_y2_ix) x(x_y1_y2_ix:x_y2_y3_ix) x(x_y2_y3_ix:x_y3_y5_ix) flip(x(1:x_y3_y5_ix))], [y1(1:x_y1_y2_ix) y2(x_y1_y2_ix:x_y2_y3_ix) y3(x_y2_y3_ix:x_y3_y5_ix) flip(y5(1:x_y3_y5_ix))], 'r')

% x=linspace(1, 5);
% y1=x.^2;
% y2=x.*(-1./4)+(9./2);
% y3=x.*0;
% isct_y1_y2 = @(x) x.^2 - (x.*(-1./4)+(9./2));
% isct_y2_y3 = @(x) x.*(-1./4)+(9./2);
% isct_y1_y3 = @(x) x.^2;
% x_y1_y2 = fzero(isct_y1_y2, 1);
% x_y2_y3 = fzero(isct_y2_y3, 1);
% x_y1_y2_ix = find(x <= x_y1_y2, 1, 'last');
% x_y2_y3_ix = find(x <= x_y2_y3, 1, 'last');
% %Fill in the Shape
% %Figure
% plot(x,y1, x,y2, x,y3); grid on;
% axis([-1,20,-1,5])
% patch([x(1:x_y1_y2_ix) x(x_y1_y2_ix:x_y2_y3_ix) flip(x(1:x_y2_y3_ix))], [y1(1:x_y1_y2_ix) y2(x_y1_y2_ix:x_y2_y3_ix) flip(y3(1:x_y2_y3_ix))], 'r')




% % lam =1.55 % micron
% % d = 0.500 % pitch in micron
% % dc = 0.6 %duty cycle
% % n_cl = 1.44
% % eps_cl=sqrt(n_cl)
% % eps_s = 1.44^2 
% % eps_r=2.6^2
% % 
% % 
% % cond1 = (Neff ~= lam/d)
% % 
% % % Neff > sqrt(eps_s)
% % cond2 = (Neff > sqrt(eps_s))
% % leakage_factor =((eps_r - eps_cl)*sin(pi*dc))^2;
% % cond3 = (Neff > n_cl) % close to n_cl
% % x=1:0.5:5;
% % y=1:0.5:5;
% 
% %[x,y]=meshgrid(1:0.5:5)
% syms x y
% % ezplot(x^2==y^4)
% n_air = 1;
% n_s=1.44
% f1=y-x-n_air; %y-x<n_air => f1<0
% f2=x-y-n_air; %x-y>n_air =f2>0
% f3=2*x-y-n_s; % 2*x-y>n_s =>f3>0
% f4=y-x; %y~=x => f4~=0
% f5=y-0.0001*x-n_s; %y>ns => f5>0
% %
% % patch(f1,f2,'r');
% ezplot(y-x==n_air,[1,5,1,5]) % f1<0
% hold on
% ezplot(x-y==n_air,[1,5,1,5]) % f2>0
% ezplot(2*x-y==n_s,[1,5,1,5]) % f3>0
% ezplot(x-y==0,[1,5,1,5])     % f4~=0
% ezplot(y-0.0001*x==n_s,[1,5,1,5]) % f5>0
% plot(f1)
% H1=area(f1);
% hold on
% idx=x>3&x<4;
% H=area(x(idx),y1(idx));
% set(H(1),'FaceColor',[1 0.5 0]);
% 
% 
% % x=0:pi/50:2*pi;
% % y1=x.^2;
% % H1=area(x,y1);
% % hold on
% % idx=x>3&x<4;
% % H=area(x(idx),y1(idx));
% % set(H(1),'FaceColor',[1 0.5 0]);