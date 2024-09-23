clc
clear
global imin imax jmin jmax nx ny sigma0 r
load('D2Data_Multy');%FVM

%绘制流线图
% figure(1)
subplot(121)
x=linspace (0 ,Lx, nx+1);
y=linspace (0 ,Ly, ny+1);
uu=u(4:imax+1,4:imax+1)';
vv=v(4:imax+1,4:imax+1)';
quiver(x,y,uu,vv);
startx = linspace(0,1,50);
starty = ones(size(startx));
streamline(x,y,uu,vv,startx,starty);
streamslice(x,y,uu,vv);
title("Re=100 streamline");
axis equal;
set(gca,'xtick',[],'ytick',[]) ; 

%绘制压力分布
% figure(2)
subplot(122)
plot(p(imin+nx/2,jmin:jmax));


normU=norm(u(imin:imax,jmin:jmax),inf);
normV=norm(v(imin:imax,jmin:jmax),inf);
sigmaCalcu=(max(p(imin+nx/2,jmin:jmax))-min(p(imin+nx/2,jmin:jmax)))*r;
error=(sigma0-sigmaCalcu)/sigma0;


fprintf(" U-∞范数=%s \n V-∞范数=%s \n Laplace Law误差=%s",num2str(normU),num2str(normV),strcat( num2str(error*100),"%"));