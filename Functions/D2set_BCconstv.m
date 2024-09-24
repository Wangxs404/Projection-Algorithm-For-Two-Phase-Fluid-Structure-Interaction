function [bu_constv] = D2set_BCconstv(av)
global velocity_BCtype
global v_bottom  v_top  v_left v_right
global D2Loc
%b-below
switch velocity_BCtype
    case 1
        bc_s=2*v_bottom *av.s.* D2Loc.s ;
        bc_n=2*v_top    *av.n.* D2Loc.n;
        bc_w=0  *av.w.* D2Loc.w;
        bc_e=0  *av.e.* D2Loc.e;
        bu_constv=bc_s+bc_w+bc_e+bc_n;
    case 2
        bc_s=2*v_bottom *av.s.* D2Loc.s ;
        bc_n=2*v_top    *av.n.* D2Loc.n;
        bc_w=2*v_left   *av.w.* D2Loc.w;
        bc_e=2*v_right  *av.e.* D2Loc.e;
        bu_constv=bc_s+bc_w+bc_e+bc_n;
    case 3
        bc_s=2*v_bottom *av.s.* D2Loc.s ;
        bc_n=2*v_top    *av.n.* D2Loc.n;
        bc_w=0  *av.w.* D2Loc.w;
        bc_e=0  *av.e.* D2Loc.e;
        bu_constv=bc_s+bc_w+bc_e+bc_n;
    case 4
        bc_s=2*v_bottom *av.s.* D2Loc.s ;
        bc_n=2*v_top    *av.n.* D2Loc.n;
        bc_w=0*v_left   *av.w.* D2Loc.w;
        bc_e=0*v_right  *av.e.* D2Loc.e;
        bu_constv=bc_s+bc_w+bc_e+bc_n;
    case 5
        bc_s=0*v_bottom *av.s.* D2Loc.s ;
        bc_n=0*v_top    *av.n.* D2Loc.n;
        bc_w=0*v_left   *av.w.* D2Loc.w;
        bc_e=0*v_right  *av.e.* D2Loc.e;
        bu_constv=bc_s+bc_w+bc_e+bc_n;
end
end
