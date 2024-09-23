function [u_star_corrected,v_star_corrected]=CorrectIBM(u_star,v_star)  
global dt A_IB EuNeighbor_Index dx DiracMatrix Arc_Length rho_Light rho_Heavy IB_Coord Coord_x Coord_y imax jmax
% IBM implicit

%计算由邻近Eu点速度构成的u*和v*向量
u_star_EuNeighbor = zeros(length(EuNeighbor_Index.X),1);
v_star_EuNeighbor = zeros(length(EuNeighbor_Index.X),1);

for i=1:length(EuNeighbor_Index.X)
    u_star_EuNeighbor(i)=u_star(EuNeighbor_Index.X(i),EuNeighbor_Index.Y(i));
    v_star_EuNeighbor(i)=v_star(EuNeighbor_Index.X(i),EuNeighbor_Index.Y(i));    
end

%构造bu和bv向量, 默认为静态结构物，即 U-IB = 0
RHS_IBM_U = 0 - dx^2 * DiracMatrix * u_star_EuNeighbor;
RHS_IBM_V = 0 - dx^2 * DiracMatrix * v_star_EuNeighbor;

% 求解线性方程组得到各拉格朗日点的加速度au和av
Accelerate_X = A_IB \ RHS_IBM_U;
Accelerate_Y = A_IB \ RHS_IBM_V;


% 将拉格朗日点的加速度，插值回临近欧拉点得到修正速度
% 遍历每各邻近Eu点，所有La点向Eu点计算分布速度并求和
u_correct = zeros(imax+3,jmax+3);
v_correct = zeros(imax+3,jmax+3);

for ii = 1:length(EuNeighbor_Index.X)

    i= EuNeighbor_Index.X(ii);
    j= EuNeighbor_Index.Y(ii);

    uEu_AuxVector = zeros(length(IB_Coord),1); % 定义辅助向量，遍历计算所有La点的Dirac分布速度并求和
    vEu_AuxVector = zeros(length(IB_Coord),1);

        for m = 1:length(IB_Coord)
            uEu_AuxVector(m) = Accelerate_X(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i),Coord_y(j)]);
            vEu_AuxVector(m) = Accelerate_Y(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i),Coord_y(j)]);
        end

    u_correct(i,j) = rho_Light*Arc_Length*dt* sum(uEu_AuxVector); % 采用rho_light ???
    v_correct(i,j) = rho_Light*Arc_Length*dt* sum(vEu_AuxVector);


end


u_star_corrected = u_star + u_correct;
v_star_corrected = v_star + v_correct;

end