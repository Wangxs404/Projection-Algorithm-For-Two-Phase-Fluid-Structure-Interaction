function [u_star_corrected, v_star_corrected] = CorrectIBM(u_star, v_star)  
% This function corrects the intermediate velocity field (u_star, v_star)
% using the Immersed Boundary Method (IBM) to account for the presence of
% solid structures within the fluid domain.

global dt A_IB EuNeighbor_Index dx DiracMatrix Arc_Length rho_Light rho_Heavy IB_Coord Coord_x Coord_y imax jmax

%==========================IBM Implicit Correction==========================

% Calculate the u* and v* vectors composed of the velocities of neighboring Eulerian points
u_star_EuNeighbor = zeros(length(EuNeighbor_Index.X), 1);
v_star_EuNeighbor = zeros(length(EuNeighbor_Index.X), 1);

for i = 1:length(EuNeighbor_Index.X)
    u_star_EuNeighbor(i) = u_star(EuNeighbor_Index.X(i), EuNeighbor_Index.Y(i));
    v_star_EuNeighbor(i) = v_star(EuNeighbor_Index.X(i), EuNeighbor_Index.Y(i));    
end

% Construct the RHS vectors for the IBM correction, assuming static structures (U-IB = 0)
RHS_IBM_U = 0 - dx^2 * DiracMatrix * u_star_EuNeighbor;
RHS_IBM_V = 0 - dx^2 * DiracMatrix * v_star_EuNeighbor;

% Solve the linear system to obtain the accelerations at Lagrangian points (au and av)
Accelerate_X = A_IB \ RHS_IBM_U;
Accelerate_Y = A_IB \ RHS_IBM_V;

% Interpolate the Lagrangian point accelerations back to the neighboring Eulerian points to obtain corrected velocities
u_correct = zeros(imax+3, jmax+3);
v_correct = zeros(imax+3, jmax+3);

for ii = 1:length(EuNeighbor_Index.X)
    i = EuNeighbor_Index.X(ii);
    j = EuNeighbor_Index.Y(ii);

    uEu_AuxVector = zeros(length(IB_Coord), 1); % Auxiliary vector to sum up Dirac distributed velocities from all Lagrangian points
    vEu_AuxVector = zeros(length(IB_Coord), 1);

    for m = 1:length(IB_Coord)
        uEu_AuxVector(m) = Accelerate_X(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i), Coord_y(j)]);
        vEu_AuxVector(m) = Accelerate_Y(m) * DiracInterpolation(IB_Coord{m}, [Coord_x(i), Coord_y(j)]);
    end

    u_correct(i, j) = rho_Light * Arc_Length * dt * sum(uEu_AuxVector); % Using rho_light ???
    v_correct(i, j) = rho_Light * Arc_Length * dt * sum(vEu_AuxVector);
end

% Apply the corrections to the intermediate velocities
u_star_corrected = u_star + u_correct;
v_star_corrected = v_star + v_correct;

end