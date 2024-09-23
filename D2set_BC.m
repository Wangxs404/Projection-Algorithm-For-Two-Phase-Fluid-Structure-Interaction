function [u,v] = D2set_BC(u,v)
%���ñ߽�����
global imax imin jmax jmin velocity_BCtype
global u_bottom v_bottom u_top  v_top  u_left  v_left  u_right v_right
%1:NS slip EW Solid_NoSlip; 2:Solid Bc ; 
%3:NS slip EW Open_NoSlip ; 4:NS slip EW Period
%%
u_top=0;    v_top=0;
u_bottom=0; v_bottom=0;
u_left=0;   v_left=0;
u_right=0;  v_right=0;

%%
switch velocity_BCtype
    case 1  
        % NS
        % U ������ �̱ڣ�k,Const��=(-1,2*u_bc)
        % V ������ �̱� (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);
        
        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        % EW
        % U ������   �̱�  (-1,2*u_bc)
        % V ���޻��� �̱�  (1,0)
        u(imin-1,:)=2*(u_left-u(imin,:))+u(imin,:);
        u(imin-2,:)=4*(u_left-u(imin,:))+u(imin,:);
        u(imin-3,:)=6*(u_left-u(imin,:))+u(imin,:);
        v(imin-1,:)=v(imin,:);
        v(imin-2,:)=v(imin,:);
        v(imin-3,:)=v(imin,:);
        
        u(imax+1,:)=2*(u_right-u(imax,:))+u(imax,:);
        u(imax+2,:)=4*(u_right-u(imax,:))+u(imax,:);
        u(imax+3,:)=6*(u_right-u(imax,:))+u(imax,:);
        v(imax+1,:)=v(imax,:);
        v(imax+2,:)=v(imax,:);
        v(imax+3,:)=v(imax,:);
        
    case 2        
        % NS
        % U ������ �̱� (-1,2*u_bc)
        % V ������ �̱� (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);
        
        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        % EW
        % U ������ �̱� (-1,2*u_bc)
        % V ������ �̱� (-1,2*v_bc)
        u(imin-1,:)=2*(u_left-u(imin,:))+u(imin,:);
        u(imin-2,:)=4*(u_left-u(imin,:))+u(imin,:);
        u(imin-3,:)=6*(u_left-u(imin,:))+u(imin,:);
        v(imin-1,:)=2*(v_left-v(imin,:))+v(imin,:);
        v(imin-2,:)=4*(v_left-v(imin,:))+v(imin,:);
        v(imin-3,:)=6*(v_left-v(imin,:))+v(imin,:);
        
        u(imax+1,:)=2*(u_right-u(imax,:))+u(imax,:);
        u(imax+2,:)=4*(u_right-u(imax,:))+u(imax,:);
        u(imax+3,:)=6*(u_right-u(imax,:))+u(imax,:);
        v(imax+1,:)=2*(v_right-v(imax,:))+v(imax,:);
        v(imax+2,:)=4*(v_right-v(imax,:))+v(imax,:);
        v(imax+3,:)=6*(v_right-v(imax,:))+v(imax,:);
        
       case 3    
        % NS
        % U ������ �̱� (-1,2*u_bc)
        % V ������ �̱� (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);
        
        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        % EW
        % U ������ ���ű߽� (1,0)
        % V ������ ���ű߽� (1,0)
        u(imin-1,:)=u(imin,:);
        u(imin-2,:)=u(imin,:);
        u(imin-3,:)=u(imin,:);
        v(imin-1,:)=v(imin,:);
        v(imin-2,:)=v(imin,:);
        v(imin-3,:)=v(imin,:);
        
        u(imax+1,:)=u(imax,:);
        u(imax+2,:)=u(imax,:);
        u(imax+3,:)=u(imax,:);
        v(imax+1,:)=v(imax,:);
        v(imax+2,:)=v(imax,:);
        v(imax+3,:)=v(imax,:);     
        
    case 4    
        % NS
        % U ������ �̱� (-1,2*u_bc)
        % V ������ �̱� (-1,2*v_bc)
        u(:,jmin-1)=2*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-1)=2*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-2)=4*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-2)=4*(v_bottom-v(:,jmin))+v(:,jmin);
        u(:,jmin-3)=6*(u_bottom-u(:,jmin))+u(:,jmin);
        v(:,jmin-3)=6*(v_bottom-v(:,jmin))+v(:,jmin);
        
        u(:,jmax+1)=2*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+1)=2*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+2)=4*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+2)=4*(v_top-v(:,jmax))+v(:,jmax);
        u(:,jmax+3)=6*(u_top-u(:,jmax))+u(:,jmax);
        v(:,jmax+3)=6*(v_top-v(:,jmax))+v(:,jmax);
        %EW Period
        %����Ϊ�Գ����ڱ߽磬��Au��B_Const�п���Case3�����޻��ƹ���
        % U ������ ���ű߽� (1,0)
        % V ������ ���ű߽� (1,0)
        u(imax+1,:)=u(imin+1,:);
        u(imax+2,:)=u(imin+2,:);
        u(imax+3,:)=u(imin+3,:);
        u(imin-1,:)=u(imax-1,:);
        u(imin-2,:)=u(imax-2,:);
        u(imin-3,:)=u(imax-3,:);
        v(imax+1,:)=v(imin+1,:);
        v(imax+2,:)=v(imin+2,:);
        v(imax+3,:)=v(imin+3,:);
        v(imin-1,:)=v(imax-1,:);
        v(imin-2,:)=v(imax-2,:);
        v(imin-3,:)=v(imax-3,:);
      case 5   %��ȫ���ڡ�ȫ�ԳƱ߽�
        % NS
        % U ������ ���ű߽� (1,0)
        % V ������ ���ű߽� (1,0)
        u(:,jmin-1)=u(:,jmin);
        v(:,jmin-1)=v(:,jmin);
        u(:,jmin-2)=u(:,jmin);
        v(:,jmin-2)=v(:,jmin);
        u(:,jmin-3)=u(:,jmin);
        v(:,jmin-3)=v(:,jmin);
        
        u(:,jmax+1)=u(:,jmax);
        v(:,jmax+1)=v(:,jmax);
        u(:,jmax+2)=u(:,jmax);
        v(:,jmax+2)=v(:,jmax);
        u(:,jmax+3)=u(:,jmax);
        v(:,jmax+3)=v(:,jmax);
        % EW
        % U ������ ���ű߽� (1,0)
        % V ������ ���ű߽� (1,0)
        u(imin-1,:)=u(imin,:);
        u(imin-2,:)=u(imin,:);
        u(imin-3,:)=u(imin,:);
        v(imin-1,:)=v(imin,:);
        v(imin-2,:)=v(imin,:);
        v(imin-3,:)=v(imin,:);
        
        u(imax+1,:)=u(imax,:);
        u(imax+2,:)=u(imax,:);
        u(imax+3,:)=u(imax,:);
        v(imax+1,:)=v(imax,:);
        v(imax+2,:)=v(imax,:);
        v(imax+3,:)=v(imax,:);     
        

        
end


end

