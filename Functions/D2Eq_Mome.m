function [u, v, uL, vL, p, pL, ppie] = D2Eq_Mome(u, v, uL, vL, uf, vf, ufL, vfL, p, mu, rho, rhoL, rhoLL, Mccx, Mccy, psi, phi, sigma, Fe, Fm)
% This function summarizes the three steps of the projection method
% Global variables are declared in subfunctions
global dt D2GradX D2GradY

% Step 1: FVM calculation of predicted velocity field [U_star, V_star]
[u_star, v_star] = D2Eq_Mome_S1(u, v, uL, vL, uf, vf, ufL, vfL, p, mu, rho, rhoL, rhoLL, Mccx, Mccy, psi, phi, sigma, Fe, Fm);
uL = u; vL = v;  % Store current velocities for BDF2 scheme

% Step 2: Apply Immersed Boundary Method correction
[u_star, v_star] = CorrectIBM(u_star, v_star);  

% FVM solution of Poisson equation to calculate ppie (pressure correction)
[u_star, v_star] = D2set_BC(u_star, v_star);  % Apply boundary conditions to predicted velocities
[p] = D2set_BCNeu(p);  % Apply Neumann boundary conditions to pressure

% Step 3: Solve for pressure correction using Rhie-Chow interpolation
[ppie] = D2RhieChow3(u_star, v_star, rho, mu, p);
[ppie] = D2set_BCNeu(ppie);  % Apply Neumann boundary conditions to pressure correction

pL = p;  % Store current pressure for next iteration
p = p + ppie;  % Update pressure for next iteration

% Step 4: Calculate U_n+1 (final velocity field)
[ppie] = D2set_BCNeu(ppie);  % Reapply Neumann boundary conditions to pressure correction
u = u_star - (2*dt/3) ./ rho .* D2GradX(ppie);  % Correct x-velocity
v = v_star - (2*dt/3) ./ rho .* D2GradY(ppie);  % Correct y-velocity
end