function [ u,v,uL,vL,p,pL,ppie] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm)
%%  此函数汇总投影法三步，仅指示输入输出即可，所用全局变量均在子函数内声明
global dt D2GradX D2GradY  
global A_IB EuNeighbor_Index dx DiracMatrix Arc_Length rho_Light IB_Coord Coord_x Coord_y imax jmax
%%
%Step1：FVM计算预测步速度场[U_star，V_star]===============================
[u_star,v_star]=D2Eq_Mome_S1(u,v,uL,vL,uf,vf,ufL,vfL,p,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);
uL=u;   vL=v;           %为引入BDF二阶，截断储存一次u，v

% %% IBM implicit
% u_star_EuNeighbor = zeros(length(EuNeighbor_Index.X),1);
% v_star_EuNeighbor = zeros(length(EuNeighbor_Index.X),1);
% 
% for i=1:length(EuNeighbor_Index.X)
%     u_star_EuNeighbor(i)=u_star(EuNeighbor_Index.X(i),EuNeighbor_Index.Y(i));
%     v_star_EuNeighbor(i)=v_star(EuNeighbor_Index.X(i),EuNeighbor_Index.Y(i));    
% end
% 
% RHS_IBM_U = 0 - dx^2 * DiracMatrix * u_star_EuNeighbor;
% RHS_IBM_V = 0 - dx^2 * DiracMatrix * v_star_EuNeighbor;
% 
% Accelerate_X = A_IB \ RHS_IBM_U;
% Accelerate_Y = A_IB \ RHS_IBM_V;
% 
% % 将拉格朗日点的加速度，插值回临近欧拉点得到修正速度
% u_correct = zeros(imax+3,jmax+3);
% v_correct = zeros(imax+3,jmax+3);
% 
% for ii = 1:length(EuNeighbor_Index.X)
% 
%     i= EuNeighbor_Index.X(ii);
%     j= EuNeighbor_Index.Y(ii);
% 
%     uEu_AuxVector = zeros(length(IB_Coord),1);
%     vEu_AuxVector = zeros(length(IB_Coord),1);
% 
%         for m = 1:length(IB_Coord)
%             uEu_AuxVector(m) = Accelerate_X(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i),Coord_y(j)]);
%             vEu_AuxVector(m) = Accelerate_Y(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i),Coord_y(j)]);
%         end
% 
%     u_correct(i,j) = rho_Light*Arc_Length*dt* sum(uEu_AuxVector);
%     v_correct(i,j) = rho_Light*Arc_Length*dt* sum(vEu_AuxVector);
% 
% 
% end
% 
% u_star = u_star + u_correct;
% v_star = v_star + v_correct;

[u_star,v_star]=CorrectIBM(u_star,v_star);
%%
%Step2：FVM解Poison方程，计算ppie=========================================
[u_star,v_star]=D2set_BC(u_star,v_star);[p]=D2set_BCNeu(p); 
[ppie] = D2RhieChow3(u_star,v_star,rho,mu,p);
[ppie] = D2set_BCNeu(ppie);%补充一次BC，否则外三层全为0，造成下次循环的P无法使用
pL=p;       %截取pL，用于计算Uf
p=p+ppie;   %基于补充压力算法时，由Ppie计算更新uv，同步更新p用于下次循环的Step1

%% 引入IBM修正

%%
% Step3：计算U_n+1========================================================
[ppie] = D2set_BCNeu(ppie);
u=u_star-(2*dt/3)./rho.*D2GradX(ppie) ; v=v_star-(2*dt/3)./rho.*D2GradY(ppie);

end

