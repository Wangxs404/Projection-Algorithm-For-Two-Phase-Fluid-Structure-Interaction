function [u_star, v_star] = D2Eq_Mome_S1(u, v, uL, vL, uf, vf, ufL, vfL, p, mu, rho, rhoL, rhoLL, Mccx, Mccy, psi, phi, sigma, Fe, Fm)
% This function calculates the intermediate velocity field (u_star, v_star)
% using the first step of the momentum equation, including various contributions
% such as advection, diffusion, pressure gradient, and external forces.

global imax imin jmax jmin gx gy nx ny Swith_SM 

%%
%==========================Initialization and Preparation===================================
u_star = zeros(imax+3, jmax+3);
v_star = zeros(imax+3, jmax+3);
uR = 2 * u - uL;
vR = 2 * v - vL;

if Swith_SM == 1 
    mX = Mccx;
    mY = Mccy;
else
    mX = rhoL .* u;
    mY = rhoL .* v;
end

[mX] = D2set_BCNeu(mX);
[mY] = D2set_BCNeu(mY);  % Momentum boundary conditions temporarily set as Neumann, not equivalent to velocity boundary conditions

%=============Calculate Coefficient Matrices=========
% Coefficient matrix A
[Au, au] = D2Matrix_Au(rho, mu);
[Av, av] = D2Matrix_Av(rho, mu);

%-----------------------
[Dtx, Dty] = D2Diver_Dtrans(uR, vR, mu); % Diffusion term 2
[bu] = D2Calcu_b(u, uL, mX, gx, @D2Gradx_Matrix, Dtx, rho, rhoL, rhoLL, uf, vf, ufL, vfL, psi, phi, p, @D2Grady_Matrix, sigma, Fe.X, Fm.X);
[bv] = D2Calcu_b(v, vL, mY, gy, @D2Grady_Matrix, Dty, rho, rhoL, rhoLL, uf, vf, ufL, vfL, psi, phi, p, @D2Gradx_Matrix, sigma, Fe.Y, Fm.Y);

%================================Calculate bu====================================
bu_final = bu - D2set_BCconstu(au);  % Constant part for ghost cells, introducing bfinal
Lu = ichol(Au);      % In some cases, no preconditioning is faster
[uv, ~] = pcg(Au, bu_final, 1e-6, 20, Lu, Lu'); % pcg is generally better than cgs
u_star(imin:imax, jmin:jmax) = reshape(uv, nx, ny);  % Assign to u_star

%================================Calculate bv====================================
bv_final = bv - D2set_BCconstv(av);
Lv = ichol(Av);
[vv, ~] = pcg(Av, bv_final, 1e-6, 20, Lv, Lv');
v_star(imin:imax, jmin:jmax) = reshape(vv, nx, ny);  % Assign to v_star

end