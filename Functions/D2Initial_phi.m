function [phi,phiL]=D2Initial_phi(initial_type)
global epsilon imax imin jmax jmin Coord_x Coord_y Ly Lx nx ny Dim r
%1.Bubble_Rising;
phi=zeros(imax+3,jmax+3);  %和u/v初始化为同样维度，方便定位,[-1 为重]

for j=jmin: jmax
    for i=imin: imax
        switch initial_type % Light is "1",Heavy is "-1"
            case 1          %Bubble_Rising
                X0=0.5*Dim;Y0=0.5*Dim;r=0.25*Dim;
                phi(i,j)  = -tanh((sqrt((Coord_x(i)-X0)^2 + (Coord_y(j)-Y0)^2)-r)/(sqrt(2)*epsilon));
        end
    end
end
phi=D2set_BCNeu(phi);   %由于未定义外层坐标，所以给一次边界条件 %可兼容润湿性，相当于初始无润湿，此后立刻润湿
phiL=phi ;
end

