function [phin, phiLn, psi, M_change] = D2Eq_PhaseLiuTVD3(uf, vf, phi)
global M epsilon dt How_Many_Phase

% Solve the Cahn-Hilliard equation discretely
% Calculate the chemical potential
psi = phi.^3 - phi - epsilon^2 * D2La_Oper(phi) + 0.0000001 * D2La_Oper(phi);
[psi] = D2set_BCNeu(psi);  % Apply Neumann boundary conditions to psi

% Calculate the mobility
M_change = M * abs((1 + phi) .* (1 - phi));

% Update phase fraction
phiLn = phi;  % Store the previous phi

switch How_Many_Phase
    case 1
        phin = phi;  % Single-phase flow: no change in phi
    case 2
        % Third-order TVD Runge-Kutta scheme for two-phase flow
        phiA = phi + dt * CHRHS(uf, vf, phi, psi, M_change);
        phiB = 0.75 * phi + 0.25 * phiA + 0.25 * dt * CHRHS(uf, vf, phiA, psi, M_change);
        phiC = (1/3) * phi + (2/3) * phiB + (2/3) * dt * CHRHS(uf, vf, phiB, psi, M_change);
        phin = phiC;
end
end

function [RHS] = CHRHS(uf, vf, phi, psi, M_change)
global D2adv D2GradX D2GradY

% Calculate diffusive flux
Difx = M_change .* D2GradX(psi);
Dify = M_change .* D2GradY(psi);

% Calculate Laplacian of chemical potential (M*∇²ψ)
MLa_Psi = D2GradX(Difx) + D2GradY(Dify);

% Calculate advective flux
Fphi = D2adv(uf, vf, uf, vf, phi, phi);  % Using lower-order time format, but with second-order time advancement

% Combine advective and diffusive terms
RHS = -Fphi + MLa_Psi;
end