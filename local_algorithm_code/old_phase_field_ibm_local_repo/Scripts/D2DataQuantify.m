% 本子程序实现数据后处理
% 方法是：在计算完成后，循环导入计算中保存的'DataNum.mat'文件，并进行相应的数据处理
%         一、避免了在计算中，每一步循环都要进行数据处理 二、便于后期增删改
%   说明：注意循环变量‘mat_Freq : mat_Freq : istep_max’
%         新增定量比数据，可不进行矩阵大小预设
clc
clear
load('D2Data_Initial.mat')
for LoadFreq = mat_Freq : mat_Freq : istep_max
    filename=strcat('D2Data',num2str(LoadFreq),'.mat');
    load(filename);
    switch initial_type           % 定量对比的数据处理
        case 1   % 计算Rising Bubble时自动开启
            phiV=(1+phi)/2  ;   LocaY=repmat(Coord_y(jmin:jmax)',1,nx)';
            Deno=sum(sum(phiV(imin:imax,jmin:jmax)));  %denominator
            NumeV=sum(sum(v(imin:imax,jmin:jmax).*phiV(imin:imax,jmin:jmax)));% Numerator V
            NumeY=sum(sum(LocaY.*phiV(imin:imax,jmin:jmax)));
            Vc(LoadFreq/mat_Freq,1)=NumeV/Deno;%#ok<SAGROW>
            Yc(LoadFreq/mat_Freq,1)=NumeY/Deno-dy/2;%#ok<SAGROW>
            
        case 4   % 计算Rayleigh Taylor时自动开启
            for j=jmin:jmax
                if phi(imin,j)*phi(imin,j+1) < 0
                    RT_wall(istep/mat_Freq)=Coord_y(j);%#ok<SAGROW>
                end
                if phi(imin+nx/2,j)*phi(imin+nx/2,j+1) < 0
                    RT_midd(istep/mat_Freq)=Coord_y(j);%#ok<SAGROW>
                end
            end
            
        case 6   % 计算剪切单液滴变形
            [col,row]=find(phi<0.55 & phi>0.45);%返回界面下标
            Lcol=(col-3)/nx*Lx;Lrow=(row-3)/ny*Ly;%将下标还原为位置坐标
            Length=zeros(size(Lcol));
            n=0;
            for j=1:size(Lcol)
                n=n+1;
                Length(n)=norm(   [1,1]   -     [Lcol(n),Lrow(n)]  );
            end
            Lmax=max(Length);Lmin=min(Length);
            Deform(istep/mat_Freq)=(Lmax-Lmin)/(Lmax+Lmin);%#ok<SAGROW>
            %体心压力梯度测试
            GPx=Grady_Matrix(p);  GPy=Grady_Matrix(p);
            DeltaPx(istep/mat_Freq)=GPx(imin+nx/2,jmin+ny/2);%#ok<SAGROW>
            DeltaPy(istep/mat_Freq)=GPy(imin+nx/2,jmin+ny/2);%#ok<SAGROW>
    end
    
end
% Time=linspace(0 , dt*istep_max , istep_max/mat_Freq);

% 由于定量化数据从第一个非零时间步开始取样，故还需在相关向量前填充一个初始向量
switch initial_type          
        case 1
          Time=dt*(mat_Freq : mat_Freq : istep_max)';
          Time=[0;Time];
          Vc=[0;Vc];%0时刻速度为0
          Yc=[0.5;Yc];%0时刻质心坐标为0.5
end

