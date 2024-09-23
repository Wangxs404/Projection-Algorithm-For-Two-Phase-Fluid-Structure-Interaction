function [IB_Coord,Euler_Index,Euler_Coord, Euler_Dirac, EuNeighbor_Index, Num_Lagrangian, DiracMatrix] = D2IBM_initial_implicit(Eu_radius, Arc_Length)
global imax imin jmax jmin Coord_x Coord_y dx

%% === IBM 拉格朗日点 ===

xo1 = 0.0; % 圆心X坐标
xo2 = 1.0;
yo = 1; % 圆心Y坐标
r_solid = 0.4;  % 半径

theta_seperated = Arc_Length / r_solid;
Num_Lagrangian = ceil(2*pi/theta_seperated);

fprintf("Number of Lagrangian Point is: %s \n\n", num2str(Num_Lagrangian));

% theta_seperated=2*pi/Num_Lag;
% arc = r_solid*theta_seperated;

IB_Coord=cell(Num_Lagrangian,1);

Vertor_IB = zeros(Num_Lagrangian,2);

for i = 1:Num_Lagrangian
    LocLag=zeros(1,2);
    if(cos(theta_seperated*(i - 1)) > 0)
        LocLag(1,1) = xo1 + r_solid*cos(theta_seperated*(i-1)) ;
        LocLag(1,2) = yo + r_solid*sin(theta_seperated*(i-1)) ;
    else
        LocLag(1,1) = xo2 + r_solid*cos(theta_seperated*(i-1)) ;
        LocLag(1,2) = yo + r_solid*sin(theta_seperated*(i-1)) ;
    end
   
    IB_Coord{i} = LocLag;

    Vertor_IB(i,1) = LocLag(1,1);
    Vertor_IB(i,2) = LocLag(1,2);

end

%% === IBM 确定邻域欧拉点 ===

Euler_Index = cell(Num_Lagrangian,1); % 表示每个La点的邻近Eu点的i、j索引

Euler_Coord = cell(2,1); % 表示每个La点的邻近Eu点的x、y坐标, 默认2行1列，后续自增

plot(Vertor_IB(:,1),Vertor_IB(:,2), 'b', 'LineWidth', 1.8);
hold on

for m=1:Num_Lagrangian
    CountEu = 0;
    IndexEuler = zeros(1,2);
    LocEuler = zeros(1,2);
    for j=jmin: jmax
        for i=imin: imax
            if sqrt((Coord_x(i)- IB_Coord{m}(1))^2 + (Coord_y(j)-IB_Coord{m}(2))^2) < Eu_radius
                CountEu = CountEu + 1;
                % 将欧拉点的索引存入元组
                IndexEuler(CountEu,1) = i;
                IndexEuler(CountEu,2) = j;
       

                % 将欧拉点的坐标存入元组
                LocEuler(CountEu,1) = Coord_x(i);
                LocEuler(CountEu,2) = Coord_y(j);
            

            end
        end
    end
    
    Euler_Index{m} = IndexEuler;
    Euler_Coord{m} = LocEuler;
    
    % fprintf("第%s个拉格朗日点，有%s个欧拉点 \n",num2str(m),num2str(CountEu));
    % fprintf("第%s个拉格朗日点，有%s个欧拉Coord \n",num2str(m),num2str(length(Euler_Coord{m})));
    % fprintf("第%s个拉格朗日点，有%s个欧拉Index \n\n",num2str(m),num2str(length(Euler_Index{m})));
    
    % 绘制每个拉格朗日点的邻近欧拉点
    scatter(Euler_Coord{m}(:,1),Euler_Coord{m}(:,2),'x');
    [x,y] = PltCicle(IB_Coord{m},Eu_radius);
    plot(x, y, 'r', 'LineWidth', 0.5);
end

%%  计算每个lagrangian点邻域Euler点的Dirac函数向量
Euler_Dirac = cell(Num_Lagrangian,1);

for m=1:Num_Lagrangian
    Dirac_Vector = zeros(length(Euler_Coord{m}),1);
    for n = 1:length(Euler_Coord{m})

        Delta_h = DiracInterpolation(Euler_Coord{m}(n,:), IB_Coord{m});
        % fprintf("第%s个拉格朗日点的第%s个欧拉点的Delta_h = %s \n\n",num2str(m),num2str(n),num2str(Delta_h));

        Dirac_Vector(n,1) = Delta_h;
        
    end
    Euler_Dirac{m} = Dirac_Vector;

end

%% 划分La点的影响域，在FIB向FEu分配时减少不必要的计算
EuNeighbor_Index.X = zeros(2,1);
EuNeighbor_Index.Y = zeros(2,1);
EuNeighbor_Coord = zeros(1,2);
EuNeighbor_cout = 0;
for j=jmin: jmax
    for i=imin: imax
        if ((sqrt((Coord_x(i)- xo1)^2 + (Coord_y(j)-yo)^2) < (r_solid + 5*dx)) && ...
           (sqrt((Coord_x(i)- xo1)^2 + (Coord_y(j)-yo)^2) > (r_solid - 5*dx))) || ...
            ((sqrt((Coord_x(i)- xo2)^2 + (Coord_y(j)-yo)^2) < (r_solid + 5*dx)) && ...
           (sqrt((Coord_x(i)- xo2)^2 + (Coord_y(j)-yo)^2) > (r_solid - 5*dx)))
            EuNeighbor_cout = EuNeighbor_cout + 1;

            EuNeighbor_Index.X(EuNeighbor_cout) = i;
            EuNeighbor_Index.Y(EuNeighbor_cout) = j;

            EuNeighbor_Coord(EuNeighbor_cout,1) = Coord_x(i);
            EuNeighbor_Coord(EuNeighbor_cout,2) = Coord_y(j);
        end
    end
end
% scatter(EuNeighbor_Coord(:,1),EuNeighbor_Coord(:,2),'c_');

%% 构造隐式IBM的Dirac矩阵
DiracMatrix=zeros(Num_Lagrangian, EuNeighbor_cout);

for m = 1:Num_Lagrangian
    for n = 1:EuNeighbor_cout
    
    DiracMatrix(m,n) = DiracInterpolation(EuNeighbor_Coord(n,:), IB_Coord{m});

    end
end



end
