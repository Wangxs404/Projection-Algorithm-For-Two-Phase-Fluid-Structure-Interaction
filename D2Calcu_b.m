function[b]=D2Calcu_b(u,uL,m,g,D2GradFunc,Dt,rho,rhoL,rhoLL,uf,vf,ufL,vfL,psi,phi,p,D2GradFunc2,sigma,Fe,Fm,F_Eu)
global ds  imin imax jmin jmax  nx ny dt D2adv dx dy  epsilon FsTension How_Many_Phase
%=============================计算bu======================================
BuTime= ds*(2*rhoL.*u-0.5*rhoLL.*uL);    
bu_time=reshape(BuTime(imin:imax,jmin:jmax),nx*ny,1);

BuGradP= -dt*ds*D2GradFunc(p);  %压力梯度项
bu_p=reshape(BuGradP(imin:imax,jmin:jmax),nx*ny,1);

BuAdv=-dt*ds*D2adv(uf,vf,ufL,vfL,m,m); %离散层面-守恒形式 U▽·(ρU）
bu_adv=reshape(BuAdv(imin:imax,jmin:jmax),nx*ny,1);%基于WENO和面心速度计算对流项，并时空积分

bu_dif=reshape(dt*(dx*dy)*Dt,nx*ny,1);   %黏性项2

BuFg= dt*ds*g*rho;  %重力项
bu_Fg=reshape(BuFg(imin:imax,jmin:jmax),nx*ny,1);

%% 表面张力 (通过 FsTension ~="ON"，关闭表面张力-sigma置零)
if FsTension ~="ON"; sigma= 0*sigma ;end 
    
if How_Many_Phase~=3  
    Fs=(0.75*sqrt(2)*epsilon)*Calcu_Fs(sigma,phi,psi,D2GradFunc,D2GradFunc2);%两相的表面张力计算
    BuFs= dt*ds*Fs;
    bu_Fs=reshape(BuFs(imin:imax,jmin:jmax),nx*ny,1);
end

b=bu_time + bu_adv + bu_p + bu_dif + bu_Fg + bu_Fs;

% 考虑F IBM
% % BuF_Eu= dt*ds*F_Eu;
% % bu_F_Eu=reshape(BuF_Eu(imin:imax,jmin:jmax),nx*ny,1);
% % b=bu_time + bu_adv + bu_p + bu_dif + bu_Fg + bu_Fs + bu_F_Eu;

end

function [Fs]=Calcu_Fs(sigma,phi,psi,D2GradFunc,D2GradFunc2)
global epsilon
Fs1= psi.*D2GradFunc(phi).*(sigma/epsilon^2); % 法向Fs
Fs2= (D2GradFunc(phi).^2 + D2GradFunc2(phi).^2) .* D2GradFunc(sigma);% 切向Fs1
Fs3= -(D2GradFunc(sigma).*D2GradFunc(phi)+D2GradFunc2(sigma).*D2GradFunc2(phi)).*D2GradFunc(phi);%切向Fs2

Fs= Fs1+Fs2+Fs3;
end
