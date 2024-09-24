function [bu_constu] = D2set_BCconstu(au)
global velocity_BCtype
global u_bottom u_top  u_left  u_right
global D2Loc
%b-below
switch velocity_BCtype
    case 1
        bc_s=2*u_bottom *au.s.* D2Loc.s ;
        bc_n=2*u_top    *au.n.* D2Loc.n;
        bc_w=2*u_left   *au.w.* D2Loc.w;
        bc_e=2*u_right  *au.e.* D2Loc.e;
        bu_constu=bc_s+bc_w+bc_e+bc_n;
    case 2
        bc_s=2*u_bottom *au.s.* D2Loc.s ;
        bc_n=2*u_top    *au.n.* D2Loc.n;
        bc_w=2*u_left   *au.w.* D2Loc.w;
        bc_e=2*u_right  *au.e.* D2Loc.e;
        bu_constu=bc_s+bc_w+bc_e+bc_n;
    case 3
        bc_s=2*u_bottom *au.s.* D2Loc.s ;
        bc_n=2*u_top    *au.n.* D2Loc.n;
        bc_w=0*u_left   *au.w.* D2Loc.w;
        bc_e=0*u_right  *au.e.* D2Loc.e;
        bu_constu=bc_s+bc_w+bc_e+bc_n;
    case 4
        bc_s=2*u_bottom *au.s.* D2Loc.s ;
        bc_n=2*u_top    *au.n.* D2Loc.n;
        bc_w=0*u_left   *au.w.* D2Loc.w;
        bc_e=0*u_right  *au.e.* D2Loc.e;
        bu_constu=bc_s+bc_w+bc_e+bc_n;
    case 5
        bc_s=0*u_bottom *au.s.* D2Loc.s ;
        bc_n=0*u_top    *au.n.* D2Loc.n;
        bc_w=0*u_left   *au.w.* D2Loc.w;
        bc_e=0*u_right  *au.e.* D2Loc.e;
        bu_constu=bc_s+bc_w+bc_e+bc_n;
end
end
