global gx gy  M epsilon  rho_Heavy rho_Light  mu_Light  mu_Heavy 
global  sigma0 

% 1. RisingBubble Case1 

switch property  %����
    case 1
        mu_Light=1  ;  mu_Heavy=10  ;  rho_Light =100  ;  rho_Heavy=1000; % ��-1��represent ��the Heavy one��
        gx=0   ;     gy=-0.98;                 % g����Ϊ��
        M=1e-2 ;     epsilon=0.02;      sigma0=24.5;        %Ǩ���� ������ ��������ϵ��
end

