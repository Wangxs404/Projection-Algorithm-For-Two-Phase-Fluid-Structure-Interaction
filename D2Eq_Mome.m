function [ u,v,uL,vL,p,pL,ppie] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm)
%%  �˺�������ͶӰ����������ָʾ����������ɣ�����ȫ�ֱ��������Ӻ���������
global dt D2GradX D2GradY  
global A_IB EuNeighbor_Index dx DiracMatrix Arc_Length rho_Light IB_Coord Coord_x Coord_y imax jmax
%%
%Step1��FVM����Ԥ�ⲽ�ٶȳ�[U_star��V_star]===============================
[u_star,v_star]=D2Eq_Mome_S1(u,v,uL,vL,uf,vf,ufL,vfL,p,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);
uL=u;   vL=v;           %Ϊ����BDF���ף��ضϴ���һ��u��v

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
% % ���������յ�ļ��ٶȣ���ֵ���ٽ�ŷ����õ������ٶ�
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
%Step2��FVM��Poison���̣�����ppie=========================================
[u_star,v_star]=D2set_BC(u_star,v_star);[p]=D2set_BCNeu(p); 
[ppie] = D2RhieChow3(u_star,v_star,rho,mu,p);
[ppie] = D2set_BCNeu(ppie);%����һ��BC������������ȫΪ0������´�ѭ����P�޷�ʹ��
pL=p;       %��ȡpL�����ڼ���Uf
p=p+ppie;   %���ڲ���ѹ���㷨ʱ����Ppie�������uv��ͬ������p�����´�ѭ����Step1

%% ����IBM����

%%
% Step3������U_n+1========================================================
[ppie] = D2set_BCNeu(ppie);
u=u_star-(2*dt/3)./rho.*D2GradX(ppie) ; v=v_star-(2*dt/3)./rho.*D2GradY(ppie);

end

