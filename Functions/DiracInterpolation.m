function [Delta_h] = DiracInterpolation(LocEu,LocLa)
global dx
% 如果欧拉点-拉格朗日点的距离超过2倍网格分辨率，则Delta_h为零，欧拉点越靠近拉格朗日点，权重越大。

h=dx;
Delta_h = (1/h^2) * Dirac((LocEu(1,1)-LocLa(1,1))/h) * Dirac((LocEu(1,2)-LocLa(1,2))/h);

end



% function [DiracDistance] = Dirac(Dis)
% % Dirac 函数
% r = abs(Dis);
% 
% if r >= 2
%     DiracDistance = 0;
% elseif r < 2 && r >= 1
%     DiracDistance = 0.125 * (5 - 2*r + sqrt( -7 + 12 * r - 4 * r^2));
% elseif r < 1 && r >= 0
%     DiracDistance = 0.125 * (3 - 2*r + sqrt( 1 + 4 * r - 4 * r^2));
% else
%     print("DiracDistance is Nan")
% end
% 
% end