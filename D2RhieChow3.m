function [ppie]=D2RhieChow3(u_star,v_star,rho,mu,p)
global  imax imin jmax jmin dt nx ny dxi dyi dx dy 
% global sigmafX sigmafY
% 中心差分计算通量+RhioChow修正
%%
%==================系数矩阵：▽·（1/ρ*▽P'）===========
ppie=zeros(imax+3,jmax+3);
[L] = D2Matrix_Laplace(rho);
%==================密度基系数
Rhox_f=2 ./ ( 1./rho(imin-1:imax , jmin:jmax) + 1./rho(imin:imax+1 , jmin:jmax) );
Rhoy_f=2 ./ ( 1./rho(imin:imax , jmin-1:jmax) + 1./rho(imin:imax , jmin:jmax+1) );


%=================其他自导系数
cx_f=3/2*Rhox_f/dt;
cy_f=3/2*Rhoy_f/dt;

%aP aE n+1)*(n+1
RCaP= 0.5*(dy/dx)*(circshift(mu,1,1)+ 2*mu + circshift(mu,-1,1))...
    +0.5*(dx/dy)*(circshift(mu,1,2)+ 2*mu + circshift(mu,-1,2));%统一计算，分别调用
RCaP_apx=RCaP(imin-1:imax ,jmin:jmax  );
RCaP_aex=RCaP(imin :imax+1,jmin:jmax  );
RCaP_apy=RCaP(imin :imax  ,jmin-1:jmax);
RCaP_any=RCaP(imin :imax  ,jmin:jmax+1);
%==========df
dx_f= 0.5*(dx*dy)*(1./RCaP_apx + 1./RCaP_aex);
dy_f= 0.5*(dx*dy)*(1./RCaP_apy + 1./RCaP_any);
%=========df cap
dx_fCap= dx_f./(1+ cx_f .* dx_f);
dy_fCap= dy_f./(1+ cy_f .* dy_f);


%系数为 dCap 
%===================b向量  ： (1.5/dt)*▽·U*===========
[FaceUX,~] = D2Matrix_FaceMeanF(u_star) ;
[~,FaceVY] = D2Matrix_FaceMeanF(v_star) ;

FaceUX_Correct=  (0.25/dx*dx_fCap)...
    .* (-1*p(imin-1-1:imax-1 , jmin:jmax)...
    +3*p(imin-1  :imax+1-1 , jmin:jmax)...
    -3*p(imin+1-1:imax+2-1   , jmin:jmax)...
    +1*p(imin+2-1:imax+3-1 , jmin:jmax));
FaceUX_RC2=FaceUX+FaceUX_Correct;

FaceVY_Correct=  (0.25/dy*dy_fCap)...
    .* (-1*p(imin:imax ,jmin-1-1:jmax-1)...
    +3*p(imin:imax ,jmin-1  :jmax+1-1)...
    -3*p(imin:imax ,jmin+1-1:jmax+2-1  )...
    +1*p(imin:imax ,jmin+2-1:jmax+3-1));
FaceVY_RC2=FaceVY+FaceVY_Correct;
adv_Ustar= dxi*( FaceUX_RC2(2:end,:)-FaceUX_RC2(1:end-1,:))...
    +dyi*(FaceVY_RC2(:,2:end)-FaceVY_RC2(:,1:end-1))  ;

Bp1=1.5/dt * adv_Ustar; b_poison=reshape(Bp1,nx*ny,1);
%%
%===================Directly Solve==============================
pv=sparse(L)\sparse(b_poison);   % Poison 方程非正定，不能使用共轭梯度法
ppie(imin : imax , jmin : jmax)=reshape(pv,nx,ny); % 将解得的P填入网格

end

