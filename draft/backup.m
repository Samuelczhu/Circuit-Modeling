%% ELEC 4700 Assignment 4 - Circuit Modeling
%
% Student: Samuel (Wendi) Zhu 
%
% Student Number: 101088968
%
% Date: 3/31/2022

% Clear all
clearvars
clearvars -global
close all
format shorte

% Make plot pretier
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth', 2);

%% Q1
% Using a fixed bottled neck value and the simulation code for assignment
% 3, do a voltage sweep of the device from 0.1V to 10V and model the
% current-voltage characteristics.

% Global variables
global C         % constants module that holds all the constants
global x y       % arrays for the current electrons positions: 1 row, colum for current position 
global xp yp     % arrays for the previous electrons positions: 1 row, column for previous position
global vx vy     % arrays for current electrons velocities: 1 row, column for current velocity
global ax ay     % scalars for electron acceleration in x and y direction
global limits    % Limits for the plot
global boxes;  % matrix for the boxes: n rows, and 4 columns for [x y w h]
% Initalize global constants
% Electron charge
C.q_0 = 1.60217653e-19;  % C         
% Rest mass
C.m0 = 9.1093837015e-31; % KG
% Effective mass of electrons
C.mn = 0.26*C.m0;  % KG
% Boltzmann constant
C.kb = 1.38064852e-23;  %m^2 * kg * s^-2 * K^-1

% Initialize the region size 200nm X 100nm
Region.x = 200e-9;
Region.y = 100e-9; 
limits = [0 Region.x 0 Region.y];  % plot limit
% Initialize the temperature
T = 300;  % K
% Initialize the mean time between collision
Tmn = 0.2e-12;  % 0.2ps
deltaT = 6e-14; % Time interval per simulation step in second
% Define the dimension
L = Region.x * 10^9;  % Length in nm
W = Region.y * 10^9;  % Width in nm 
% boxLF = 0;
% boxWF = 0;
boxLF = 0.3;  % Fraction of the length of the box
boxWF = 0.4;  % Fraction of the width of the box
Lb = boxLF*L;  % Length of the box in nm
Wb = boxWF*W;  % Width of the box in nm
deltaXY = 0.02*L;  % Assume deltaX = deltaY in nm
% Calculate the dimension of solution matrix
nx = (L/deltaXY);
ny = (W/deltaXY);

% % Number of sweep points for the voltage
% nSweep = 30;
% vecVoltage = linspace(0.1, 10, nSweep);  % Generate the voltage vector
% vecCurrent = zeros(1, nSweep);  % vector for holding the current

% Voltage - deltaT
%  10     -  3e-14
%  5      -  6e-14

% % Loop through the different voltages
% for iVolt = 1:length(vecVoltage)
    % Define the voltages
%     voltageX0 = vecVoltage(iVolt);  % Voltage at X=0
%     voltageX1 = 0;  % Voltage at X=L
    
voltageX0 = 10;
voltageX1 = 0;

    % Step 1: Calculate the E field
    % Calculate the meshgrid
    [X,Y] = meshgrid(linspace(0,L,nx), linspace(0,W,ny));    
    % Declare the matrix for conductivity: Sigma(y,x)
    matrixSigma = ones(ny, nx);  % Dimension: ny times nx
    xIndexBox = ceil((L-Lb)/(2*deltaXY));  % Find the starting x index for the box
    LbIndexRange = ceil(Lb/deltaXY);  % Index range for the length of the box
    WbIndexRange = ceil(Wb/deltaXY);  % Index range for the width of the box
    % Assign the region for the box
    matrixSigma(1:WbIndexRange, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;
    matrixSigma(ny-WbIndexRange:ny, xIndexBox:xIndexBox+LbIndexRange) = 10^-2;    
    % Declare the matrix for voltage V(y,x)
    matrixV = zeros(ny, nx);  % Dimension: ny times nx   
    % Declare the G matrix and F vector: GV = F
    G = zeros(nx*ny, nx*ny);  
    F = zeros(nx*ny, 1);  
    % Construct the G matrix and F vector
    for ix = 1:nx
        for iy = 1:ny
            % Calculate the index
            n = mappingEq(ix, iy, ny);
            % Check for the boundary
            if ix==1 || ix==nx || iy ==1 || iy==ny
                G(n,n) = 1;
                % Boundary condition for x
                if ix == 1
                    F(n,1) = voltageX0;  % V at x = 0 
                elseif ix == nx
                    F(n,1) = voltageX1;  % and V at x = L
                elseif iy == 1
                    nyp = mappingEq(ix, iy+1, ny);  % dV/dy=0 at y=0
                    G(n,nyp) = -1;
                elseif iy == ny
                    nym = mappingEq(ix, iy-1, ny);  % dV/dy=0 at y=W
                    G(n, nym) = -1;
                end
            else
                % Calculate the sigma
                sigmaxp = (matrixSigma(iy,ix) + matrixSigma(iy,ix+1))/2;
                sigmaxm = (matrixSigma(iy,ix) + matrixSigma(iy, ix-1))/2;
                sigmayp = (matrixSigma(iy,ix) + matrixSigma(iy+1, ix))/2;
                sigmaym = (matrixSigma(iy,ix) + matrixSigma(iy-1, ix))/2;         
                % Calculate mapping index
                nxp = mappingEq(ix+1, iy, ny);  % index for V(i+1,j)
                nxm = mappingEq(ix-1, iy, ny);  % index for V(i-1,j)
                nyp = mappingEq(ix, iy+1, ny);  % index for V(i,j+1)
                nym = mappingEq(ix, iy-1, ny);  % index for V(i,j-1)    
                % Setup the G matrix
                G(n,n) = -(sigmaxp+sigmaxm+sigmayp+sigmaym)/deltaXY^2;
                G(n, nxp) = sigmaxp/deltaXY^2;
                G(n, nxm) = sigmaxm/deltaXY^2;
                G(n, nyp) = sigmayp/deltaXY^2;
                G(n, nym) = sigmaym/deltaXY^2;
            end
        end
    end
    % Solve for V from GV = F
    V = G\F;
    % Map back to the 2D region
    for iMap = 1:nx*ny
        % Calculate the index for the 2D region
        ix = ceil(iMap/ny);
        iy = mod(iMap, ny);
        if iy == 0
            iy = ny;
        end
        % Assign the value
        matrixV(iy, ix) = V(iMap);
    end
    % Solve the electric field
    [Ex, Ey] = gradient(-matrixV);
    Ex = Ex/(deltaXY * 10^-9);  % convert to V/m
    Ey = Ey/(deltaXY * 10^-9);  % convert to V/m

    % Step 2: Calculate the acceleration field
    % Initialize the number of "super" electrons
    numE = 1000;
    % Number of simulation steps
    numSim = 1000;
    % Boudary mode: specular(0) or diffusive(1)
    boundaryMode = 0;
    % Add the boxes
    numBox = AddObstacles(boxLF, boxWF, Region);
    % To find the current, the following steps are performed:
    % 1) Calculate the total area
    areaA = Region.x * Region.y;  % m^2
    areaA = areaA * 100^2;  % cm^2
    % 2) Calculate the total electrons in the area assuming electron
    % concentration is 10^15 cm-2
    totalE = 10^15 * areaA;  % total electrons
    % 3) Find the charge per "Super Electron", where "Super Electron" is the
    % particle in this simulation
    numEPerSuperE = totalE/numE;  % number of electron per super electron
    superECharge = -C.q_0 * numEPerSuperE;  % Charge per super electron
    % 4) The current can be found by counting the net number of super electrons

    % Initialize acceleration for each electron
    ax = zeros(1, numE);  % Acceleration in x
    ay = zeros(1, numE);  % Acceleration in y
    % Calculate the acceleration field: a = Force/mass = q*E/mass
    accFieldX = -C.q_0 * Ex / (C.mn);
    accFieldY = -C.q_0 * Ey / (C.mn);

    % Add the electrons
    AddElectrons_WithBox(numE, Region, T, numBox);
    % Calculate the scattering probability
    Pscat = 1-exp(-deltaT/Tmn);
    
    % TODO: delete later
    numEPlot = 5;
    numGridX = 15; % number of grid in x direction
    numGridY = 15; % number of grid in y direction
    % Initalize plot
    figure(1)
    axCol = axes;
    axCol.ColorOrder = rand(numEPlot,3); % Initalize color for each electron
    hold on

    % Super electron count for current calculation
    % Count on left side x=0. +1 flow right, -1 flow left
    countECurrent = 0;  % Hold the super electron count

    %% TODO delete later
    vecCurr = zeros(1, numSim);
    vecTime = zeros(1, numSim);
    timeIndex = 0;
    % Step 3: Loop for simulation
    for iSim = 1:numSim
        % TODO: delete later 
%         PlotPoint(numEPlot, numGridX, numGridY);

        % Store the current positions
        xp = x;
        yp = y;
        % Calculate the future positions: x = x0 + vx*t
        x = x + vx * deltaT;
        y = y + vy * deltaT;
        % Calculate the future velocity: vx = ax*t
        vx = vx + ax*deltaT;
        vy = vy + ay*deltaT;

        % Reset the super electron count
        countECurrent = 0;  

        % Loop through all the particles 
        for iE=1:numE
            % flag for invalid position
            bInvalid = false;
            % Step 1 - Check for boundary
            % Check for invalid x position
            if x(iE) <= 0
                x(iE) = Region.x; % Appear on right
                xp(iE) = x(iE); 
                bInvalid = true;
                % Update the electron count for current calculation
                countECurrent = countECurrent-1;  % -1 flow left
            elseif x(iE) >= Region.x
                x(iE) = 0; % Appear on left 
                xp(iE) = x(iE);
                bInvalid = true;
                % Update the electron count for current calculation
                countECurrent = countECurrent+1;  % +1 flow right
            end
            % Check for invalid y position
            if y(iE) <= 0
                bInvalid = true;
                y(iE) = 0;
                % Check for boundary mode
                if boundaryMode == 0  % Specular boundary
                    vy(iE) = -vy(iE);
                else  % Diffusive boundary  TODO: check diffusive implementation
                    vy(iE) = abs(sqrt(C.kb*T/C.mn).*randn());  % positive vy
                end
            elseif y(iE) >= Region.y
                y(iE) = Region.y;
                bInvalid = true;
                % Check for boundary mode
                if boundaryMode == 0  % Specular boundary
                    vy(iE) = -vy(iE);
                else  % Diffusive boundary
                    vy(iE) = -abs(sqrt(C.kb*T/C.mn).*randn());  % negative vy
                end
            end   
            % Step 2: Check for boxes
            for iBox = 1:numBox
                % Retrieve box info
                boxX1 = boxes(iBox, 1);
                boxX2 = boxes(iBox, 1)+boxes(iBox, 3);
                boxY1 = boxes(iBox, 2);
                boxY2 = boxes(iBox, 2)+boxes(iBox, 4);
                % Check if the particle is inside a box
                if (x(iE)>=boxX1 && x(iE)<=boxX2 && y(iE)>=boxY1 && y(iE) <= boxY2)
                    bInvalid = true;   %Invalid position
                    % Check for x position
                    if xp(iE) <= boxX1  % Coming from left side
                        x(iE) = boxX1;
                        % Check for boundary mode
                        if boundaryMode == 0  % Specular boundary
                            vx(iE) = -vx(iE);
                        else  % Diffusive boundary
                            vx(iE) = -abs(sqrt(C.kb*T/C.mn).*randn());  % negative vx
                        end                    
                    elseif xp(iE) >= boxX2  % Coming from right side
                        x(iE) = boxX2;                   
                        % Check for boundary mode
                        if boundaryMode == 0  % Specular boundary
                            vx(iE) = -vx(iE);
                        else  % Diffusive boundary
                            vx(iE) = abs(sqrt(C.kb*T/C.mn).*randn());  % positive vx
                        end
                    end
                    % Check for y position
                    if yp(iE) <= boxY1  % Coming from bottom
                        y(iE) = boxY1;               
                        % Check for boundary mode
                        if boundaryMode == 0  % Specular boundary
                            vy(iE) = -vy(iE);
                        else  % Diffusive boundary
                            vy(iE) = -abs(sqrt(C.kb*T/C.mn).*randn());  % negative vy
                        end
                    elseif yp(iE) >= boxY2 % Coming from top
                        y(iE) = boxY2;                   
                        % Check for boundary mode
                        if boundaryMode == 0  % Specular boundary
                            vy(iE) = -vy(iE);
                        else  % Diffusive boundary
                            vy(iE) = abs(sqrt(C.kb*T/C.mn).*randn());  % positive vy
                        end
                    end
                    % Break the loop for box
                    break;
                end
            end            
            % Step 3: Check for scattering
            if ~bInvalid && Pscat > rand()
                % Rethermalize  TODO: Check rethermalize process is correct
                vx(iE) = sqrt(C.kb*T/C.mn).*randn();
                vy(iE) = sqrt(C.kb*T/C.mn).*randn();
            end    
            % Step 4: Find acceleration
            % Find the corresponding index for the acceleration field
            indexX = ceil(x(iE)/(deltaXY*10^-9));
            indexY = ceil(y(iE)/(deltaXY*10^-9));
            % Check for invalid index
            if indexX <= 0
                indexX = 1;
            end
            if indexY <= 0
                indexY = 1;
            end
            % Assign the acceleration of the electron
            ax(iE) = accFieldX(indexX);
            ay(iE) = accFieldY(indexY);
        end
        %% TODO: delete later
        % Pause some time
%         pause(0.01);
        timeIndex = timeIndex + deltaT;
        vecCurr(iSim) = superECharge*countECurrent/deltaT;
        vecTime(iSim) = timeIndex;
    end

    % Calculate the current
%     vecCurrent(iVolt) = superECharge*countECurrent/deltaT;

figure(2)
plot(vecTime, vecCurr, "-b.")
    
% end

% % Plot the current versus voltage characteristics
% figure(1)
% plot(vecVoltage, vecCurrent, "-b.")
% title("Current - Voltage Characteristics")
% xlabel("Voltage (V)")
% ylabel("Current (A)")
% % Do the linear fit for determining resistance value of the device
% polyR = polyfit(vecVoltage, vecCurrent, 1);
% grid on



















%% Helper functions
% The following functions are helper functions used in the main code

% Helper function for mapping index
% @param iRow = i index for the row
%        jRow = j index for the column
%        ny   = size of the y
function [n] = mappingEq(iRow, jCol, ny)
    n = jCol + (iRow - 1) * ny;
end  % End mappingEq


% Helper function to add the obstacles
% @ param  boxLF = length of the box in fraction of region.x
%          boxWF = width of the box in fraction of region.y
%          region = region.x and region.y
function [numBox] = AddObstacles(boxLF, boxWF, region)
global boxes  % Matrix for holding the boxes
% Find the x, y, w, h for the bottom box
xbb = region.x/2 - region.x*boxLF/2;
ybb = 0;
wbb = region.x*boxLF;
hbb = region.y * boxWF;
% Find the x, y, w, h for the upper box
xub = region.x/2 - region.x * boxLF/2;
yub = region.y * (1-boxWF);
wub = region.x * boxLF;
hub = region.y * boxWF;
% Create the boxes
boxes = [xbb ybb wbb hbb;
    xub yub wub hub];
% Return number of boxes
numBox = height(boxes);
end  % End AddObstacles


% This function add a bunch of electrons in a given region randomly for Q3
% @param numE = number of electrons
%        region = region for the electrons
%        T = temperature in Kelvin
%        numBox = number of boxes
function AddElectrons_WithBox(numE, region, T, numBox)
global C  % Constants
global x y % arrays for current electron positions
global xp yp % arrays for previous electron positions
global vx vy % arrays for current electron velocities
global boxes  % Matrix for the boxes position

% Create the arrays for electrons locations
x = rand(1, numE) * region.x;
y = rand(1, numE) * region.y;

% Loop through the electrons to make sure that no electrons inside obstacles
for iE = 1:numE
    % Flag to indicate whether inside box
    insideBox = true;
    while (insideBox)
        insideBox = false;
        % Loop through the boxes
        for iBox = 1:numBox
            % Check for invalid electrons position
            if (x(iE)>boxes(iBox, 1) && x(iE)<(boxes(iBox, 1)+boxes(iBox, 3)) ...
                    && y(iE)>boxes(iBox, 2) && y(iE) < (boxes(iBox, 2)+boxes(iBox, 4)))
                insideBox = true;
                break;
            end
        end
        if (insideBox)
            % Regenerate position
            x(iE) = rand() * region.x;
            y(iE) = rand() * region.y;
        end
    end
end
% Create the arrays for previous electron positions
xp = x; 
yp = y;
% Create helper arrays for velocity distrubution
vx = sqrt(C.kb*T/C.mn).*randn(1, numE);
vy = sqrt(C.kb*T/C.mn).*randn(1, numE);
end  % End AddElectrons_WithBox



% TODO: delete later
% Helper function to plot one point for electron position
% @param  numEPlot = number of electrons to be plotted
%         numGridX = number of grid on the x axis
%         numGridY = number of grid on the y axis
function PlotPoint(numEPlot, numGridX, numGridY)
global x y xp yp limits

% plot the electron positions
plot([xp(1:numEPlot);x(1:numEPlot)], [yp(1:numEPlot);y(1:numEPlot)])

% Adjust the axis limits
axis(limits)
% Set grid
set(gca,'xtick',linspace(0, limits(2), numGridX));
set(gca, 'ytick',linspace(0, limits(4), numGridY));
grid on

% Add title and labels
title("Electron Modeling");
xlabel("Length (m)")
ylabel("Width (m)")
end  % End PlotPoint

