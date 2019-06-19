% % input: dimension  => output: frequency, wavelength and modes
% 
% clear all;
% close all;
% a=0.40e-6;
% b=0.22e-6;
% neff=2.86;
% c=300e6;
% lambda_c=0;
% i=0;
% for n=0:2
%     for m=0:2
%         if ((m || n)==1)
%             i=i+1;
%             lambda_c=(2/sqrt((m/a)^2+(n/b)^2)); %micron
%             f_c=c*sqrt(((m*pi)/a)^2+((n*pi)/b)^2)/(2*pi*neff);
%             lambda_c_prim=(c/(neff*f_c));
%             lambda_0=lambda_c*neff;
%             f_c*lambda_c*neff
%             f_id(i)=f_c;
%             lamda_id(i)=1e6*lambda_0;
%             m_id(i)=m;
%             n_id(i)=n;
%             sprintf('m=%d n=%d ==> lambda_c =%2.2e f_c = %2.2e lambda_c_prim =%2.2e lambda_0 =%2.2f ',m,n,lambda_c,f_c, lambda_c_prim,lambda_0)
%         end        
%     end
% end
% f_id
% lamda_id
% [f_id_sorted,I_sorted]=sort(f_id)
% lamda_id_sorted=lamda_id(I_sorted)
% m_id_sorted=m_id(I_sorted)
% n_id_sorted=n_id(I_sorted)



% %%%%%%%%%%%%%%%%%
% input: dimension  => output: frequency, wavelength and modes 
%sweeping width (a)

clear all;
close all;
fileId=fopen('designSpecParam.txt','w');
%a=0.40e-6;
b=0.22e-6;
neff=2.86;
c=300e6;
lambda_c=0;
set=0;
for neff=2.5:0.1:3.1
for a=0.20e-6:0.100e-6:0.60e-6
    i=0;
    set=set+1;
    for n=0:2
        for m=0:2
            if ((m || n)==1)
                i=i+1;
                lambda_c=(2/sqrt((m/a)^2+(n/b)^2)); %micron
                f_c=c*sqrt(((m*pi)/a)^2+((n*pi)/b)^2)/(2*pi*neff);
                lambda_c_prim=(c/(neff*f_c));
                lambda_0=lambda_c*neff;
                f_c*lambda_c*neff;
                f_id{set}(i)=f_c;
                lamda_id{set}(i)=1e6*lambda_0;
                m_id{set}(i)=m;
                n_id{set}(i)=n;
                %sprintf('m=%d n=%d ==> lambda_c =%2.2e f_c = %2.2e lambda_c_prim =%2.2e lambda_0 =%2.2f ',m,n,lambda_c,f_c, lambda_c_prim,lambda_0)
            end        
        end
    end
    f_id{set};
    lamda_id{set};
    [f_id_sorted{set},I_sorted{set}]=sort(f_id{set});
    lamda_id_sorted{set}=lamda_id{set}(I_sorted{set});
    m_id_sorted{set}=m_id{set}(I_sorted{set});
    n_id_sorted{set}=n_id{set}(I_sorted{set});
    lambda_min = lamda_id_sorted{set}(2);
    lambda_max = lamda_id_sorted{set}(1);
    m_min =m_id_sorted{set}(1);
    m_max =m_id_sorted{set}(2);
    n_min =n_id_sorted{set}(1);
    n_max =n_id_sorted{set}(2);
    a_set(set)=a;
    neff_set(set)=neff;
    sprintf('set=%d a=%2.2f , neff =%2.2f => lambda window %2.2f - %2.2f, Mode_mn_min=%d%d Mode_mn_max=%d%d', set, a*1e6, neff, lambda_min, lambda_max, m_min,n_min,m_max,n_max)
    fprintf(fileId,'set=%d a=%2.2f  , neff =%2.2f => lambda window %2.2f - %2.2f, Mode_mn_min=%d%d Mode_mn_max=%d%d\n', set, a*1e6, neff, lambda_min, lambda_max, m_min,n_min,m_max,n_max);
end
end

