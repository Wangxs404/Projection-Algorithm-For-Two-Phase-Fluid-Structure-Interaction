function [b] = D2Calcu_b(u, uL, m, g, D2GradFunc, Dt, rho, rhoL, rhoLL, uf, vf, ufL, vfL, psi, phi, p, D2GradFunc2, sigma, Fe, Fm, F_Eu)
% This function calculates the source term 'b' for the momentum equation
% using various contributions including time derivative, pressure gradient,
% advection, diffusion, gravity, and surface tension.

global ds imin imax jmin jmax nx ny dt D2adv dx dy epsilon FsTension How_Many_Phase

%=============================Calculate bu======================================
BuTime = ds * (2 * rhoL .* u - 0.5 * rhoLL .* uL);    
bu_time = reshape(BuTime(imin:imax, jmin:jmax), nx*ny, 1);

BuGradP = -dt * ds * D2GradFunc(p);  % Pressure gradient term
bu_p = reshape(BuGradP(imin:imax, jmin:jmax), nx*ny, 1);

BuAdv = -dt * ds * D2adv(uf, vf, ufL, vfL, m, m); % Discrete-level conservative form U▽·(ρU)
bu_adv = reshape(BuAdv(imin:imax, jmin:jmax), nx*ny, 1); % Convection term based on WENO and face-centered velocities, integrated in time and space

bu_dif = reshape(dt * (dx * dy) * Dt, nx*ny, 1);   % Diffusion term

BuFg = dt * ds * g * rho;  % Gravity term
bu_Fg = reshape(BuFg(imin:imax, jmin:jmax), nx*ny, 1);

%% Surface Tension (Turned off if FsTension is not "ON", sigma set to zero)
if FsTension ~= "ON"; sigma = 0 * sigma; end 
    
if How_Many_Phase ~= 3  
    Fs = (0.75 * sqrt(2) * epsilon) * Calcu_Fs(sigma, phi, psi, D2GradFunc, D2GradFunc2); % Surface tension calculation for two phases
    BuFs = dt * ds * Fs;
    bu_Fs = reshape(BuFs(imin:imax, jmin:jmax), nx*ny, 1);
end

b = bu_time + bu_adv + bu_p + bu_dif + bu_Fg + bu_Fs;

end

function [Fs] = Calcu_Fs(sigma, phi, psi, D2GradFunc, D2GradFunc2)
% This function calculates the surface tension force 'Fs'

global epsilon

Fs1 = psi .* D2GradFunc(phi) .* (sigma / epsilon^2); % Normal Fs
Fs2 = (D2GradFunc(phi).^2 + D2GradFunc2(phi).^2) .* D2GradFunc(sigma); % Tangential Fs1
Fs3 = -(D2GradFunc(sigma) .* D2GradFunc(phi) + D2GradFunc2(sigma) .* D2GradFunc2(phi)) .* D2GradFunc(phi); % Tangential Fs2

Fs = Fs1 + Fs2 + Fs3;

end