function [u_IB, v_IB] = IB_Velocity(u, v, Euler_Dirac,Euler_Index,Num_Lagrangian)
global dx

u_IB = zeros(Num_Lagrangian,1);
v_IB = zeros(Num_Lagrangian,1);

for m=1:Num_Lagrangian

    u_factor = zeros(length(Euler_Index{m}),1);
    v_factor = zeros(length(Euler_Index{m}),1);

    for k=1:length(Euler_Index{m})

        % Euler点速度向量
        i = Euler_Index{m}(k,1);
        j = Euler_Index{m}(k,2);
        
        u_factor(k) = u(i,j);
        v_factor(k) = v(i,j);

    end

    % 点积 
    uIB_AuxVector = dx^2 * u_factor .* Euler_Dirac{m};
    vIB_AuxVector = dx^2 * v_factor .* Euler_Dirac{m};

    % 求和
    u_IB(m) = sum(uIB_AuxVector);
    v_IB(m) = sum(vIB_AuxVector);

end


end

