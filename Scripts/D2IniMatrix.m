global  sigma0 %sigmafX sigmafY
%初始化矩阵存储
u=zeros(imax+3,jmax+3) ; v=zeros(imax+3,jmax+3);
uL=zeros(imax+3,jmax+3); vL=zeros(imax+3,jmax+3);
p=zeros(imax+3,jmax+3) ; ppie=zeros(imax+3,jmax+3) ;
% Matoolu=zeros(imax+3,jmax+3) ; Matoolv=zeros(imax+3,jmax+3);     %容器，更新面心速度使用
% [u,v]=D2set_BC(u,v) ; 

% u=zeros(imax+3,jmax+3)+1;

uf=D2fXmean(u);vf=D2fYmean(v);

ufL=uf ;  vfL=vf; %插值面心速度
psi=zeros(imax+3,jmax+3);


T=zeros(imax+3,jmax+3) ;
%适用于全局为均匀sigma时，便于使用RhieChow Fs修正，此处等同于Nuem边界
%sigma非均匀时在相关子函数中更新sigma(温度方程|表活方程)
if How_Many_Phase==3; sigma0=0 ;end
sigma=sigma0*ones(imax+3,jmax+3); 

% sigmafX=D2fXmean(sigma);sigmafY=D2fYmean(sigma);

%surfascrtan
% cs=zeros(imax+3,jmax+3);

PhiE=zeros(imax+3,jmax+3) ;%便于在关闭电场方程时，export2可以正常运行
ChargeQ=zeros(imax+3,jmax+3) ;
Fe.X=zeros(imax+3,jmax+3) ;
Fe.Y=zeros(imax+3,jmax+3) ;

PhiM= zeros(nx+6,ny+6);%初始化%便于在关闭磁场方程时，export2可以正常运行
Fm.X= zeros(imax+3,jmax+3) ;
Fm.Y= zeros(imax+3,jmax+3) ;

F_Eu.X= zeros(imax+3,jmax+3) ;%初始化%便于在关闭IBM方程时，export2可以正常运行
F_Eu.Y= zeros(imax+3,jmax+3) ;