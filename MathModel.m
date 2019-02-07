%% Foreword
% This is the main script for battery pack and motor sizing information for
% the PROVE Project 2 Endurance EV. The script generates sizing information
% and graphs for different cell types, drag, Crr, and velocity values. Any
% questions should be forwarded to Ryan Sadeghi (ryan.sadeghi@gmailcom)
% hey
% Latest Revision
% February 7th, 2019 

%% To Do

% Add acceleration loss 
    % Include number of times accelerating/stopping
% Weight distribution affect
% Total battery capacity vs time 
    % Voltage drop effect
    % Resistance of pack
    % Initial voltage effects 
% Redo empty weight analysis
    
%% Considerations

% CFD info
    % Density effect on drag
    % Velocity effect on drag
    
% Crr change 
    % how tire pressure changes over the run
    % deformation

% Powertrain efficiency 
    % Motor efficiency
    % Motor controller efficiency
    % Resistance of pack 
    
%% Clear the Deck

clc; clear all; close all;

%% Input Margins
velocity_margin = 0;
crr_margin = 0;
empty_mass_margin = 0;
drag_margin = 0;
downforce_margin = 0;
range_margin = -50; 

%% Final Car Assumed Specs

% Ambient conditions
g = 9.81;                                       % Gravitational acceleration[m/s^2]
rho = 1.225;                                    % Sea level air density [kg/m^3]

% CFD Inputs
d = 190+drag_margin;                                        % Drag force from CFD as of 11/15/18 [N]
downforce = (1100+downforce_margin)/9.81;                         % [N to kg]

% Vehicle specs
crr = .0068;                                    % Coefficient of rolling resistance value of BMW i3 Ecopia tire (model??) from Japanese manufacturer via Bridgestone
empty_mass = (1500+empty_mass_margin)*0.453592 + downforce;         % Empty mass of car based on estimations found on Confluence > Motor Sizing [kg]
v = 65+velocity_margin;                                         % Velocity [mph]

range = 1100 + range_margin;

%% Cell Specs
% Pulled from spec sheets of most viable commerically available batteries

% INR21700 M50 21700 
cell_2170.cost = 6.5;                                                   % Estimate of cost per cell commercially [$]
cell_2170.V = 3.63;                                                     % Nominal cell voltage [V]
cell_2170.energy = 18.2;                                                % Nominal cell energy [Wh]
cell_2170.capacity = cell_2170.energy/cell_2170.V;                      % Nominal cell capacity [Ah]
cell_2170.mass = 68;                                                    % Cell mass [g]
cell_2170.diameter = 21.1;                                              % Approximate diameter [mm] 
cell_2170.height = 70.15;                                               % Approximate height [mm] 
cell_2170.volume = (pi*cell_2170.diameter^2)/4 * cell_2170.height;      % Approximate volume [mm^3]
cell_2170.en_density = cell_2170.energy/cell_2170.mass;                 % Energy density [Wh/g]
cell_2170.energy_cost = cell_2170.energy/cell_2170.cost;                % Energy density cost [Wh/dollar] 

%% Cell Packing Efficiency
% Approximate maximum packing efficiency if packed in a staggered pattern (90%)
 
cell_2170.volume = cell_2170.volume/.9; 

%% Concise Power Calculations

calc = batterysize(cell_2170,empty_mass,d,v,range,crr);


fprintf('2170 Baseline Cell Outputs \n\n')
fprintf('Total weight of cells: %1.3f lbs \n', calc.num_cells(end)*cell_2170.mass/1000*2.20462)
fprintf('Total weight of car: %1.3f lbs\n', (calc.m(end)-downforce)*2.20462)
fprintf('Total number of cells: %1.3f \n', calc.num_cells(end))
fprintf('Capacity of the battery pack: %1.3f kWh\n', calc.kwh(end))
fprintf('Total cell volume: %1.3f ft^3 \n', calc.num_cells(end)*cell_2170.volume*3.53147e-8)
fprintf('Power: %1.3f kW \n', calc.P(end))
fprintf('Total cost of cells: %1.3f dollars \n', calc.num_cells(end)*cell_2170.cost)

%% Functions


function calc = batterysize(cell,empty_mass,drag,vel,range,crr)
g = 9.81;

vel = vel*.44704;
t = range/(vel/.44704); 

i = 1; 

% First iteration calculates empty car values                                                       
calc.m(i) = empty_mass;                              % Sets mass to empty car weight [kg]

calc.num_cells(i) = 0;

calc.Fd(i) = drag + crr*g*calc.m(i)*cosd(2.9);                 % force required to overcome drag [N]

calc.P(i) = (calc.Fd(i)*vel)/1000;                   % power (kW)
calc.P(i) = calc.P(i)/.85;                           % factor in powertrain losses

calc.kwh(i) = calc.P(i)*t;                           % capacity
calc.need(i) = calc.kwh(i);                          % how much extra capacity needed

while calc.need(i) > eps
    i = i+1;
    % add a cell to subtract energy needed
    calc.m(i) = calc.m(1) + cell.mass/1000*(i-1); % add weight of cell to total mass (kg)
    
    calc.num_cells(i) = calc.num_cells(i-1)+1; % add to cell counter
    
    calc.Fd(i) = drag + crr*g*calc.m(i)*cosd(2.9); % force required to overcome drag (Newtons)
    
    calc.P(i) = (calc.Fd(i)*vel)/1000; % power (kW)
    calc.P(i) = calc.P(i)/.85; % factor in powertrain losses
    
    calc.kwh(i) = calc.P(i)*t; % capacity
    calc.need(i) = calc.kwh(i) - cell.energy/1000*(i-1); % how much extra capacity needed    
end
end
