%=================Print 无量纲群
DDL=1 ;  DDU=0.01; DDrho=rho_Light  ;  DDmiu=mu_Light;         %特征长度和特征速度
fprintf('   CFL= %s  Cn= %s  \n   Pe = %s    Re= %s    Ca= %s \n\n',...
            num2str(DDU*dt/dx),       num2str(epsilon/DDL),   ...
            num2str(DDU*epsilon/M), num2str(DDrho*DDU*DDL/DDmiu), num2str(DDmiu*DDU/sigma0));%Ca数有待查证
  
