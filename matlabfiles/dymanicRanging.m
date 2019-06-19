clear all
d=0.55
L=100;
neff_a=1.9
lam_a_1=1.15
lam_a_2=1.4
n_cl=1
theta_a_1=asind(((neff_a - (lam_a_1)/d))/n_cl)
theta_a_2=asind(((neff_a - (lam_a_2)/d))/n_cl)
%%%%%%%%%%%%%%
% fwhm= (lam^2/(2*pi*L*Neff))*1000 
fwhm_a_1= (lam_a_1^2/(2*pi*L*neff_a))*1000
fwhm_a_2= (lam_a_2^2/(2*pi*L*neff_a))*1000

fwhm_b_2 = 0.6*fwhm_a_2

% fwhm_b_1 / fwhm_a_1 = (lam_b_1/lam_a_1)^2 * (neff_a / neff_b)
temp1 =(1/lam_a_2)^2 * (fwhm_a_2/fwhm_b_2) *neff_a 

% neff_b=temp1*lam_b_1^2
% neff_b = (lam_b_1/lam_a_1)^2 * (fwhm_a_1/fwhm_b_1) *neff_a 
% 
% 
% % Now we go to active ranging
% % theta_b_1 = theta_a_1
% % ((neff_a - (lam_a_1)/d))/n_cl = ((neff_b - (lam_b_1)/d))/n_cl
% 
lam_b_1 =1.15
neff_b = neff_a + (lam_b_1 - lam_a_1)/d
lam_b_2 =(neff_b - neff_a)*d +lam_a_2 

fwhm_a_1= (lam_a_1^2/(2*pi*L*neff_a))*1000
fwhm_a_2= (lam_a_2^2/(2*pi*L*neff_a))*1000

fwhm_b_1= (lam_b_1^2/(2*pi*L*neff_b))*1000
fwhm_b_2= (lam_b_2^2/(2*pi*L*neff_b))*1000

% neff_b = neff_a + (lam_b_2 - lam_a_2)/d
% syms neff_b lam_b_1 lam_b_2
% 
% eq1 = neff_b-temp1*lam_b_2^2 ==0;
% eq2 = neff_b -neff_a - (lam_b_1 - lam_a_1)/d ==0;
% eq3 = neff_b - neff_a - (lam_b_2 - lam_a_2)/d ==0;
% sol = solve([eq1, eq2, eq3], [neff_b, lam_b_1, lam_b_2]);
% Xsol = vpa(sol.neff_b)
% Ysol = vpa(sol.lam_b_1)
% Zsol = vpa(sol.lam_b_2)