%% Two-Phase Flow Solver with Finite Volume Method
% This program solves incompressible two-phase flow using a finite volume method
% with WENO and a phase field model, based on collocated grid Rhie-Chow interpolation
% and a projection algorithm.
% Copyright by Wangxs

clear; close all; warning off;

%% Global Variables
global imax imin jmax jmin nx ny Coord_x Coord_y dx dy dxi dyi dt Lx Ly Swith_Mcc
global Convection_Scheme refere_Pressure velocity_BCtype theta FsTension
global Dim initial_type dat_Freq Dynamic_Draw mat_Freq Res_Freq RhieChow_or_Not How_Many_Phase Datpath
global A_IB EuNeighbor_Index DiracMatrix Arc_Length rho_Light IB_Coord

% Import necessary functions and scripts
addpath('./Functions/');
addpath('./Scripts/');

%% Grid & Time Step Setup
Dim = 1;                % Characteristic length scale
Lx = 1 * Dim;           % Width of computational domain
Ly = 2 * Dim;           % Height of computational domain
nx = 64;                % Number of grid cells in x-direction
ny = 128;                % Number of grid cells in y-direction
t_total = 3;            % Total simulation time
dt = 1e-4;              % Time step size
istep_max = ceil(t_total / dt);  % Maximum number of time steps

% Default settings
theta = 90;           % Default wetting angle
FsTension = "Off";    % Surface tension control switch

% Multi-physics selection ("ON" to enable)
Wetting = "Off";
Thermal = "Off";
Surfactant = "Off";
Electronic = "Off";
Magnetic = "Off";
ImmerseBoundaryCondition = "Off";

% Debug and output settings
Res_Freq = 50;        % Frequency of residual display (every 50 time steps)
dat_Freq = 1000;      % Frequency of data output (every 1000 time steps)
mat_Freq = 10000;     % Frequency of matrix output (every 10000 time steps)
Dynamic_Draw = 'ON';  % Enable dynamic drawing during simulation

% Mesh generation
[imin, imax, jmin, jmax, Coord_x, Coord_y, dx, dy, dxi, dyi] = D2Mesh(Lx, Ly);

%% Convection Scheme & Interface Selection
How_Many_Phase = 2;       % Handles single, two, and three-phase flows
Convection_Scheme = 1;    % Convection scheme: 1:WENO, 2:Quick, 3:Centre
refere_Pressure = 4;      % Pressure reference point: 1/2/3/4/5 (bottom-left counterclockwise, center)
velocity_BCtype = 1;      % Check 'D2set_BC'
initial_type = 1;         % See 'D2Initial_phi'
RhieChow_or_Not = 3;      % See 'D2Eq_Mome_S2_RhieChow'
Swith_Mcc = 0;
property = 1;             % See 'D2Property'

% Run property setup script
run D2Property;           % Specify material, interface

%% Data Storage Setup
name = "RisingBubble_IBM_3";  % Create folder in parent directory to store data
Matpath = strcat("Va1-", name, "-Mat");
Datpath = strcat("Va1-", name, "-Data");
mkdir('./Output/', strcat("Va1-", name, "-Data"));
mkdir('./Output/', strcat("Va1-", name, "-Mat"));

%% Preprocessing
% Initialize phase interface
[phi, phiL] = D2Initial_phi(initial_type);

% IBM Lagrangian Points Initialization
Arc_Length = 1.2 * dx;        % Referenced from Zhu Yang's PhD thesis
Eu_radius = 1.8 * sqrt(2) * dx;  % Euler neighborhood radius, each Lagrangian point has about 19~21 Euler points
[IB_Coord, Euler_Index, Euler_Coord, Euler_Dirac, EuNeighbor_Index, Num_Lagrangian, DiracMatrix] = D2IBM_initial_implicit(Eu_radius, Arc_Length);
A_IB = rho_Light * Arc_Length * dt * dx^2 * (DiracMatrix * DiracMatrix');

% Initialize istep for continuation calculations before "SpecialSet"
istep = 0;

% Run setup scripts
run D2LocEleBc;     % Set up wall boundary element position vectors
run D2HandleFunc;   % Initialize operator function handles
run D2IniMatrix;    % Initialize matrix sizes
run D2SpecialSet;   % Handle continuation/case-specific settings
run D2DimLess;      % Set up dimensionless parameters for plotting and output
run D2PostIni;      % Specify dimensionless quantities and output initial data

%% Main Computation Loop
% The time-stepping and solution procedures starts here. Subsequent code would implement

% Close any open figures from preprocessing
close 

% Initialize residuals and timing variables
Residual = zeros(istep_max, 6);
R_limit = 100;
time_start = clock;
A_timeStart = datestr(now, 'HH:MM mmm.dd');

% Main simulation loop
while (istep < istep_max)
    % Set Boundary Conditions
    [u, v] = D2set_BC(u, v);
    [p] = D2set_BCNeu(p);
    phi = D2set_BCNeu(phi);  % Default no wetting
    if Wetting == "ON"
        phi = D2set_BCWet(phi, phiL, sigma);  % If wetting is enabled, modify ¦Õ_BC
    end
 
    % Phase Field Equations
    phiLL = phiL;
    [phin, phiLn, psi, M_change] = D2Eq_PhaseLiuTVD3(uf, vf, phi);
    phi_R = norm(phin - phi) / (nx * ny);  % Calculate residual
    phiL = phiLn;
    phi = phin;
    [rho, rhoL, rhoLL, mu, phi, phiL, phiLL] = D2Updata_RhoMiu(phi, phiL, phiLL, sigma, Wetting);  % Update ¦Ñ, ¦Ì, and constrain ¦Õ [-1,1]
    [Mccx, Mccy] = D2Calcu_Mcc(rhoL, u, v, M_change, psi);  % Correct mass flux Mcc

    % Navier-Stokes Equations
    [un, vn, uLn, vLn, pn, pL, ppien] = D2Eq_Mome(u, v, uL, vL, uf, vf, ufL, vfL, p, mu, rho, rhoL, rhoLL, Mccx, Mccy, psi, phi, sigma, Fe, Fm);

    % Calculate residuals and update variables
    u_R = norm(un - u) / (nx * ny);
    v_R = norm(vn - v) / (nx * ny);
    p_R = norm(pn - p) / (nx * ny);
    u = un; v = vn; uL = uLn; vL = vLn; p = pn; ppie = ppien;
    [uf, vf, ufL, vfL] = D2Uf_inte(u, v, uL, vL);  % Interpolate to calculate Uf

    % End of single iteration

    % Display residuals dynamically
    istep = istep + 1;
    Residual(istep, 1) = u_R;
    Residual(istep, 2) = v_R;
    Residual(istep, 3) = phi_R;
    
    if mod(istep, Res_Freq) == 0
        fprintf('\n Iteration %s / %s\n  phi_R = %s\n  u_R = %s\n  v_R = %s\n', ...
            num2str(istep), num2str(istep_max), num2str(phi_R), num2str(u_R), num2str(v_R));      
        if Thermal == "ON"
            Residual(istep, 4) = T_R;
            fprintf('  T_R = %s\n', num2str(T_R));
        end
        if Surfactant == "ON"
            Residual(istep, 5) = cs_R;
            fprintf('  cs_R = %s\n', num2str(cs_R));
        end
    elseif (u_R > R_limit) || (v_R > R_limit)
        fprintf('Residuals too large, not converging\n');
        break;
    end

    % Output and Draw
    if mod(istep, dat_Freq) == 0 || istep == 1
        % Specify output variables and ranges
        var1 = phi(imin:imax, jmin:jmax);
        var2=psi(imin:imax,jmin:jmax);
        var3=u(imin:imax,jmin:jmax);
        var4=v(imin:imax,jmin:jmax);
        var5=p(imin:imax,jmin:jmax);
        var6=T(imin:imax,jmin:jmax);
        var7=cs(imin:imax,jmin:jmax);
        var8=cpsi(imin:imax,jmin:jmax);
        var9=sigma(imin:imax,jmin:jmax);
        var10=PhiE(imin:imax,jmin:jmax);
        var11=ChargeQ(imin:imax,jmin:jmax);
        var12=PhiM(imin:imax,jmin:jmax);
        varName = 'phi psi u v p T cs cpsi sigma PhiE ChargeQ PhiM\n';
        D2export(istep, xx(:), yy(:), NGX, NGY, varName, var1(:), var2(:), var3(:), var4(:), var5(:), var6(:), var7(:), var8(:), var9(:), var10(:), var11(:), var12(:));
    end

    % Dynamic drawing
    switch Dynamic_Draw
        case 'ON'
            if mod(istep, 20) == 0
                contourf(Coord_x(4:1+nx), Coord_y(4:1+ny), phi(4:1+nx, 4:1+ny)', 3);
                title('Moving phi'); axis equal;
                drawnow
            end
    end

    % Save data to Workspace
    if mod(istep, mat_Freq) == 0
        filename = strcat('D2Data', num2str(istep));
        MatData = char(strcat('../', Matpath, '/', filename));
        save(MatData);
    end

    % Predict total runtime at 100 iterations
    if istep == 100
        time_end100 = clock;
        runtime100 = etime(time_end100, time_start) / 60;
        PredictTime = runtime100 * istep_max / 100;
        fprintf(' PredictTime = %s \n', num2str(PredictTime));
    end
end

% End of simulation, record timing information
time_end = clock;
A_timeEnd = datestr(now, 'HH:MM mmm.dd');
A_runtime_s = etime(time_end, time_start);
A_runtime_min = etime(time_end, time_start) / 60;
save(char(strcat('../', Matpath, '/', 'D2Data_Multy')));
fprintf(' EndStep= %s \n ',num2str(istep));





