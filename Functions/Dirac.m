function [DiracDistance] = Dirac(Dis)
% Dirac 函数
r = abs(Dis);

if r >= 2
    DiracDistance = 0;
elseif r < 2 && r >= 1
    DiracDistance = 0.125 * (5 - 2*r + sqrt( -7 + 12 * r - 4 * r^2));
elseif r < 1 && r >= 0
    DiracDistance = 0.125 * (3 - 2*r + sqrt( 1 + 4 * r - 4 * r^2));
else
    print("DiracDistance is Nan")
end

end