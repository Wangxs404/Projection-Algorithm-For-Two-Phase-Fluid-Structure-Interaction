function [F_IB] = IB_Force(U_IB,V_IB,IB_Coord)
global dt dx

% F_IB.X = zeros(Num_Lagrangian,1);
% F_IB.Y = zeros(Num_Lagrangian,1);

k = 1; % 直接力法，默认k=1

% 设置默认参数，调试状态下为标量，后续应改为矢量
U_IBtarget = 0;
V_IBtarget = 0;
rho_scalar = 1;

F_IB.X = k * (2*rho_scalar/dt) * (U_IB - U_IBtarget);
F_IB.Y = k * (2*rho_scalar/dt) * (V_IB - V_IBtarget);


% for m= 1:length(IB_Coord)
% quiver(IB_Coord{m}(1),IB_Coord{m}(2), F_IB.X(m), F_IB.Y(m), 1, 'LineWidth', 2); 
% hold on
% end

end