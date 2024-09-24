function [rho, rhoL, rhoLL, mu, phi, phiL, phiLL] = D2Updata_RhoMiu(phi, phiL, phiLL, sigma, Wetting)
global Reconstrct rho_Heavy rho_Light mu_Light mu_Heavy

% Update density and viscosity based on the phase field (phi)
% phi, phiL, phiLL represent the phase field at current, previous, and two time steps ago

% Constrain phi, phiL, and phiLL to the range [-1, 1]
range1 = (phi > 1); range2 = (phi < -1);
phi(range1) = 1; phi(range2) = -1;

range3 = (phiL > 1); range4 = (phiL < -1);
phiL(range3) = 1; phiL(range4) = -1;

range5 = (phiLL > 1); range6 = (phiLL < -1);
phiLL(range5) = 1; phiLL(range6) = -1;

% Apply Neumann boundary conditions (default: no wetting)
[phi] = D2set_BCNeu(phi);
[phiL] = D2set_BCNeu(phiL);
[phiLL] = D2set_BCNeu(phiLL);

% Apply wetting boundary conditions if enabled
if Wetting == "ON"
    [phi] = D2set_BCWet(phi, phiL, sigma);
    [phiL] = D2set_BCWet(phiL, phiLL, sigma);
    [phiLL] = D2set_BCWet(phiLL, phiLL, sigma);
end

% Reconstruct density and viscosity fields
rho = Reconstrct(rho_Light, rho_Heavy, phi);
rhoL = Reconstrct(rho_Light, rho_Heavy, phiL);
rhoLL = Reconstrct(rho_Light, rho_Heavy, phiLL);
mu = Reconstrct(mu_Light, mu_Heavy, phi);  % Only need mu at current time step

% Apply Neumann boundary conditions to viscosity (needed for Rhie-Chow interpolation)
[mu] = D2set_BCNeu(mu);
end