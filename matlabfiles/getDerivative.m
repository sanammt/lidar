clear all
close all
%syms neff
% lam=1.55;
% d=0.55;
% ncl=1.44;
% neff = [1.8 2.2 2.4 2.6 2.9 3]
% % neff = 3.2
% for neff= [1.8 2.2 2.4 2.6 2.9 3]
%     x=sin((neff-(lam/d))/ncl)
%     y=diff(x)
% end
%     plot(neff,y)
    
% neff=1.8:0.1:3.4 
% y=0.6944*neff-1.9571;
% x=asind(0.6944*neff-1.9571)
% plot(neff,y)

%lam1=1.35:0.05:1.65;
syms Neff
lam1=1.45
lam2=1.65
% Neff=1.8:0.1:2.75
d=0.55
n_cl=1.0
theta1 = vpa(asind((Neff-(lam1/d))/n_cl));
theta2 = vpa(asind((Neff-(lam2/d))/n_cl));
fov=theta1 - theta2;
Neff=[1.8 2.2 2.4 2.6 2.8 3.2];
z1=subs(theta1)
z2=subs(theta2)
z=subs(fov)
plot(Neff,z1);
hold on
plot(Neff,z2)
plot(Neff,z)


clear all
ncl=1.44;
% % idx=0
% idx_neff_d_set=0
idx_neff=0;
idx_t=0;
idx_neff_d_set=0;
idx_notValid = 0;
L=100; %micron
for Neff=2:0.2:2.95
%     disp('new Neff')
    idx_neff=idx_neff+1;
    idx_d=0;
    for d=0.55:0.05:0.7
%         disp('new d')
        idx_d=+idx_d+1;
        idx_neff_d_set = idx_neff_d_set+1;
        idx_lam = 0;
        for lam=1.35:0.05:1.65
            %disp('new lam-----------')
            idx_lam=idx_lam+1;
            idx_t=idx_t+1;
            Neff;
            d;
            lam;
            neff_d_set(idx_neff_d_set, 1)=Neff;
            neff_d_set(idx_neff_d_set, 2)=d;
            theta(idx_t)=asind(Neff-(lam/d))/ncl;
            theta_air(idx_t)=asind(Neff-(lam/d));
            fwhm(idx_t) = (lam^2/(2*pi*L*Neff))*1000; %nano meter
            if (~isreal(theta(idx_t)))
                idx_notValid = idx_notValid +1;
                theta(idx_t)
                notValidSet(idx_notValid,1)=Neff;
                notValidSet(idx_notValid,2)=d;
                notValidSet(idx_notValid,3)=lam;
                notValidSet(idx_notValid,4)=idx_t;
            end
            devThetatoLam(idx_t)=-1/(d*sqrt( ncl^2- (Neff-(lam/d))^2 ));
            idx_neff;
            idx_d;
            idx_lam;
            idx_t;
            idx_neff_d_set;
        end
            Neff;
            d;
            fov(idx_neff_d_set)=theta(idx_t) - theta(idx_t-idx_lam+1);
            fov(idx_neff_d_set);
            fov_air(idx_neff_d_set)=theta_air(idx_t) - theta_air(idx_t-idx_lam+1);
            fov_air(idx_neff_d_set);
            
            fwhm_worst(idx_neff_d_set) = max(fwhm(idx_t-idx_lam+1:idx_t));
            sprintf('Neff=%2.2f , d=%2.2f , fov=%2.2f, fov_air=%2.2f, fwhm_worst=%2.2f',Neff,d,fov(idx_neff_d_set), fov_air(idx_neff_d_set), fwhm_worst(idx_neff_d_set))
    end
end




% y=vpa(diff(theta1,Neff))
% Neff=[2.2 2.4 2.6];
% z=subs(y)
% plot(z)





