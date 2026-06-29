function [u,v] = D2InitialUVdb(u,v,Coord_x,Coord_y,imax,imin,jmax,jmin)
%…Ë÷√db's ÀŸ∂»
global ny nx imin imax jmin jmax
u=zeros(imax+3,jmax+3) ; v=zeros(imax+3,jmax+3);
uL=zeros(imax+3,jmax+3); vL=zeros(imax+3,jmax+3);

u_top=1;    u_bottom=-1;
uTvec=linspace(0,u_top,ny/2);uinT=repmat(uTvec,nx,1);
uBvec=linspace(u_bottom,0,ny/2);uinB=repmat(uBvec,nx,1);
u(imin:imax,jmin:jmax)=[uinB,uinT];
v(imin:imax,jmin:jmax)=0;

end

