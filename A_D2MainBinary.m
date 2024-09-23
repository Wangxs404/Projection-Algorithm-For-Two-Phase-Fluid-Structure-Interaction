%                               ˵    ��
%  ���������ͬλ����RhieChow��ֵ��ͶӰ�㷨���������������+WENO+�ೡģ����ⲻ��ѹ������
%  Copy right by Wangxs


clear ; close all ; warning off
%=================================Global Variable============================
global imax imin jmax jmin nx ny Coord_x Coord_y dx dy dxi dyi dt Lx Ly Swith_Mcc
global Convection_Scheme refere_Pressure  velocity_BCtype  theta FsTension 
global Dim initial_type dat_Freq  Dynamic_Draw mat_Freq Res_Freq RhieChow_or_Not How_Many_Phase  Datpath
global A_IB EuNeighbor_Index DiracMatrix Arc_Length rho_Light IB_Coord


%==================================��������&ʱ�䲽================================
Dim=1;
Lx=1*Dim    ;    Ly=2*Dim  ;         nx=64;    ny=128;
dt=1e-4     ;    istep_max=ceil(20/dt);  theta=90;  % Ĭ����ʪ��
FsTension="Off";    %���ƿ���Fs
Wetting="Off"; Thermal="Off";  Surfactant="Off";      % ѡ���������"ON"��
Electronic="Off"; Magnetic="Off"; ImmerseBoundaryCondition="ON";
Res_Freq=50  ;  dat_Freq=1000;  mat_Freq=10000;  Dynamic_Draw='ON'; % Debug����
[imin,imax,jmin,jmax,Coord_x,Coord_y,dx,dy,dxi,dyi] = D2Mesh(Lx,Ly);

%===============================������ʽ�����ѡ��=================
How_Many_Phase=2;                      %��˵��ࡢ���ࡢ����������
Convection_Scheme=1;                   %������ʽ  1:WENO 2:Quick 3:Centre
refere_Pressure=4;                     %ѹ�����  1/2/3/4/5��������ʱ���ĶԽ�-����
velocity_BCtype=1;  % Check��D2set_BC��   %1:NS slip EW Solid_NoSlip; 2:Solid Bc ; 3:NS slip EW Open_NoSlip ; 4:NS slip EW Period 5��All Period
initial_type=1;     % See��D2Initial_phi' %1.Bubble_Rising;2.Dam_Break;3.Zalesak's Disk;4.Rayleigh Taylor------
Swith_Mcc=2;        % See��D2Calcu_MccSu' %����ͨ�� 1:ON
RhieChow_or_Not=3;  % See��D2Eq_Mome_S2_RhieChow' %1.Centre;9. Full RhieChow-LJY;11. Full RhieChow-HZY
property=1;         % See��D2Property'    %���ϲ��� (�����Property��)
run D2Property      %===ָ�����ϡ����桢������������Ҫ��λ��IniPhi֮ǰ��

%==================================�����ݴ洢���ϼ�Ŀ¼=============
name="RisingBubble_IBM_3";        % ���ϼ�Ŀ¼�����ļ��д洢����
Matpath=strcat("Va1-",name,"-Mat");
Datpath=strcat("Va1-",name,"-Data");
mkdir('../',strcat("Va1-",name,"-Data"))  
mkdir('../',strcat("Va1-",name,"-Mat"))

%==============================ǰ����==============================
[phi,phiL]=D2Initial_phi(initial_type) ; %===������ʼ��=====

% === IBM �������յ� ��ʼ��===
Arc_Length = 1.2*dx;         % �ο�Zhu Yang�Ĳ�ʿ����
Eu_radius = 1.8*sqrt(2)*dx;  % ����Euler�뾶,ÿ��La��Լ19~21��Eu��
[IB_Coord, Euler_Index, Euler_Coord, Euler_Dirac, EuNeighbor_Index, Num_Lagrangian,DiracMatrix] = D2IBM_initial_implicit(Eu_radius, Arc_Length);
A_IB = rho_Light * Arc_Length * dt * dx^2 * (DiracMatrix * DiracMatrix');

istep=0;                                 %===Ϊ�������㣬�轫��ʼistep���� SpecialSet��ǰ
 %=====�߽�Ԫλ������========���Ӻ������=============��ʼ�������С====
 %====����/Cases���Ի�=====��ͼ�����/�����ʼ����%===ָ���������������
        run D2LocEleBc  ;   run D2HandleFunc ;  run D2IniMatrix            
        run D2SpecialSet;    run D2DimLess ;  run D2PostIni
        
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  ��ʼ�ջ�����  $$$$$$$$$$$$$$$$$$$$$$$$$$$
% close %�ر�ǰ����ͼ��

Residual=zeros(istep_max,6);R_limit=100;  time_start=clock ; A_timeStart=datestr(now,'HH:MM mmm.dd');
while (istep < istep_max)
    % ��ʼ��������  
    %=========================Set Boudary Condition==================
    [u,v]=D2set_BC(u,v)    ;    [p]=D2set_BCNeu(p);  
    phi = D2set_BCNeu(phi);      % Ĭ������ʪ��
    if Wetting=="ON" ; phi = D2set_BCWet(phi,phiL,sigma);end %��������ʪ�ԣ����޸Ħ�_BC
 
    %==============================Phase Field================
     phiLL=phiL; 
    [phin,phiLn,psi,M_change] = D2Eq_PhaseLiuTVD3(uf,vf,phi);
     phi_R=norm(phin-phi)/(nx*ny) ; phiL=phiLn; phi=phin;                             % ����в� ����ֵphi & phiL
    [rho,rhoL,rhoLL,mu,phi,phiL,phiLL]=D2Updata_RhoMiu(phi,phiL,phiLL,sigma,Wetting); % ���¦ѡ��̣�ͬʱԼ����[-1,1]
    [Mccx,Mccy]=D2Calcu_Mcc(rhoL,u,v,M_change,psi);                                   % ��������ͨ��Mcc

    %==============================MultiPhysics======================

    %==============================NS Equation=========================
    [un,vn,uLn,vLn,pn,pL,ppien] = D2Eq_Mome(u,v,uL,vL,uf,vf,ufL,vfL,p,mu,rho,rhoL,rhoLL, Mccx, Mccy,psi,phi,sigma,Fe,Fm);

    u_R = norm(un-u)/(nx*ny);    v_R=norm(vn-v)/(nx*ny);  p_R = norm(pn-p)/(nx*ny);
    u=un;   v=vn;  uL=uLn; vL=vLn;  p=pn;   ppie=ppien;
    [uf,vf,ufL,vfL]=D2Uf_inte(u,v,uL,vL);                   %���ò�ֵ��������Uf .��ѡһ

    %==================================��ɵ���ѭ��============================
    %% ��̬��ʾ�в�
    istep=istep+1 ; Residual(istep,1)=u_R ;  Residual(istep,2)=v_R ; Residual(istep,3)=phi_R ;  %����в�
    if   mod(istep,Res_Freq) ==  0
        fprintf('\n ��������  %s / %s\n  phi_R=   %s\n  u_R  =   %s\n  v_R  =   %s  \n',...
            num2str(istep),num2str(istep_max),num2str(phi_R),num2str(u_R),num2str(v_R));      
        if  Thermal=="ON"
            Residual(istep,4)=T_R ;
            fprintf('  T_R  =   %s\n',num2str(T_R));
        end
        if  Surfactant=="ON"
            Residual(istep,5)=cs_R ;
            fprintf('  cs_R =   %s\n',num2str(cs_R));
        end
    elseif (u_R>R_limit)||(v_R>R_limit)
        fprintf('�в�������� \n');
        break;
    end
    %% Output and Draw
%=====������ݵ�.dat======
if  mod(istep,dat_Freq) ==  0 || istep==1
    var1=phi(imin:imax,jmin:jmax); %1.ָ�������������Χ
    var2=psi(imin:imax,jmin:jmax);
    var3=u(imin:imax,jmin:jmax);
    var4=v(imin:imax,jmin:jmax);
    var5=p(imin:imax,jmin:jmax);
    var6=T(imin:imax,jmin:jmax);
    var7=cs(imin:imax,jmin:jmax);
    var8=cpsi(imin:imax,jmin:jmax);
    var9=sigma(imin:imax,jmin:jmax);
    var10=PhiE(imin:imax,jmin:jmax);
    var11=ChargeQ(imin:imax,jmin:jmax);
    var12=PhiM(imin:imax,jmin:jmax);
    varName ='phi psi u v p T cs cpsi sigma PhiE ChargeQ PhiM\n';  % 2.ָ�����������ʾ��;ע��ͬһ.m�ļ��к����varName�ᱻǰ��ĸ���
    D2export(istep,xx(:),yy(:),NGX,NGY,varName,var1(:),var2(:),var3(:),var4(:),var5(:),var6(:),var7(:),var8(:),var9(:),var10(:),var11(:),var12(:))%3.�������
end
%=====���ƶ�̬ͼ======
switch Dynamic_Draw
    case 'ON'
        if     mod(istep,20) ==  0
            contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),phi(4:1+nx,4:1+ny)',3) %�ȸ���ͼ
            % contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),u(4:1+nx,4:1+ny)') %�ȸ���ͼ
            % Vorticity=D2Gradx_Matrix(v)-D2Grady_Matrix(u);
            % contourf(Coord_x(4:1+nx),Coord_y(4:1+ny),Vorticity(4:1+nx,4:1+ny)',20)   %����ͼ

            title('Moving phi');axis equal;
            drawnow
        end
        % fprintf('U_IB: %s\n', mat2str(U_IB));
end

%=====�������ݵ�Workspace====
if   mod(istep,mat_Freq) ==  0   %�������ݵ�Workspace
    filename=strcat('D2Data', num2str(istep));
    MatData=char(strcat('../',Matpath,'/',filename));%��·���ͱ�����ȫ��ƴ�ӳ�char���ͣ���save����
    save (MatData);
end
if  istep==100   %�������ݵ�Workspace
    time_end100=clock;  runtime100=etime(time_end100,time_start)/60;
    PredictTime=runtime100*istep_max/100;
    fprintf(' PredictTime = %s \n',num2str(PredictTime));
end
end

time_end=clock  ; A_timeEnd=datestr(now,'HH:MM mmm.dd');
A_runtime_s=etime(time_end,time_start); A_runtime_min=etime(time_end,time_start)/60;
save(char(strcat('../',Matpath,'/','D2Data_Multy')))
% run D2DataQuantify
% fprintf(' EndStep= %s \n ',num2str(istep));





