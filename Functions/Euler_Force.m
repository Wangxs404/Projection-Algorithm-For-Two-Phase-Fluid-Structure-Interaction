function [F_Eu] = Euler_Force(F_IB, EuNeighbor_Index,IB_Coord, arc)
global imax jmax Coord_x Coord_y dx
% 反推全局Eu点的F_Eu,在流固边界附近划分影响域，提升计算效率。

F_Eu.X = zeros(imax+3,jmax+3) ;
F_Eu.Y = zeros(imax+3,jmax+3) ;

for ii = 1:length(EuNeighbor_Index.Y)

    i= EuNeighbor_Index.X(ii);
    j= EuNeighbor_Index.Y(ii);

    uEu_AuxVector = zeros(length(IB_Coord),1);
    vEu_AuxVector = zeros(length(IB_Coord),1);

        for m = 1:length(IB_Coord)
            uEu_AuxVector(m) = F_IB.X(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i),Coord_y(j)]) * dx * arc;
            vEu_AuxVector(m) = F_IB.Y(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i),Coord_y(j)]) * dx * arc;
        end

    F_Eu.X(i,j) = sum(uEu_AuxVector);
    F_Eu.Y(i,j) = sum(vEu_AuxVector);


end

end