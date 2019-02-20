%% Foreword
% This is the main script for battery pack and motor sizing information for
% the PROVE Project 2 Endurance EV. The script generates sizing information
% and graphs for different cell types, drag, Crr, and velocity values. Any
% questions should be forwarded to Ryan Sadeghi (ryan.sadeghi@gmailcom)

% Latest Revision
% February 17th, 2019 

clc; clear all; close all;

%% To Do

% Add acceleration loss 
    % Include number of times accelerating/stopping
% Weight distribution affect
% Total battery capacity vs time 
    % Voltage drop effect
    % Resistance of pack
    % Initial voltage effects 
    
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

%% Motor Information
RPM = [64 110 134 146 160 171 183 193 202 211 221 233 249 275 285 295 308 323 335 ...
    350, 366 386, 405, 425, 445, 468, 494, 516, 545, 570, 592, 627, 672, 711, ...
    766, 835, 856, 862, 870, 875, 879, 883, 887, 892, 896, 900, 902, 904, 909, ...
    913, 918, 923, 925, 928, 930, 932, 933];

T = [279.6, 278, 275.5, 273.2, 270.8, 268, 265.8, 264.2, 261.5, 258.7, 255.7, 252.6 ...
    248.2, 239.8, 235.7, 231.6, 227.2, 221.8, 216.9, 212.1, 207.3, 200.8, 195.1 ...
    189.5, 183.6, 178, 171.3, 165.6, 158.5, 152.2, 147.8, 142.9, 136, 129.7, 122.8 ...
    114.7, 107.5, 100.7, 91.2, 83.2, 77, 70.1, 64.2, 56.9, 51.4, 46, 40, 34.3, 29.2 ...
    23.9, 18.5, 11.8, 7.9, 5.5, 3.4, 2.3, .7];

% Find y-intercept
p = polyfit(RPM(1:2),T(1:2),1);
y = polyval(p,0);

RPM = [0 RPM];
T = [y T];

% Linearly interpolate between each point

i = 1;
warning off;

for i = 1:length(RPM)-1
    p(i,1:2) = polyfit(RPM(i:i+1),T(i:i+1),1);
end

%% Final Car Assumed Specs

% Ambient conditions
g = 9.81;                                       % Gravitational acceleration[m/s^2]
rho = 1.225;                                    % Sea level air density [kg/m^3]
                                       % Velocity [mph]

% CFD Inputs
d = 190+drag_margin;                                        % Drag force from CFD as of 2/9/19 [N]
downforce = (1100+downforce_margin)/9.81;                         % [N to kg]
v = 65+velocity_margin; 

% Vehicle specs
crr = .0068;                                                   % Coefficient of rolling resistance value of BMW i3 Ecopia tire (model??) from Japanese manufacturer via Bridgestone
empty_mass = (2000+empty_mass_margin)*0.453592;    % Empty mass of car based on estimations found on Confluence > Motor Sizing [kg]

range = 1050 + range_margin;

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


%% Concise Power Calculations

calc = batterysize(cell_2170,empty_mass,v,range,crr,p);

fprintf('2170 Baseline Cell Outputs \n\n')
fprintf('Total weight of cells: %1.3f lbs \n', calc.num_cells(end)*cell_2170.mass/1000*2.20462)
fprintf('Total weight of car: %1.3f lbs\n', (calc.m(end)-downforce)*2.20462)
fprintf('Total number of cells: %1.3f \n', calc.num_cells(end))
fprintf('Capacity of the battery pack: %1.3f kWh\n', calc.kwh(end))
fprintf('Total cell volume: %1.3f ft^3 \n', calc.num_cells(end)*cell_2170.volume*3.53147e-8)
fprintf('Power: %1.3f kW \n', calc.P(end))
fprintf('Total cost of cells: %1.3f dollars \n', calc.num_cells(end)*cell_2170.cost)

%% Carpet Plots

speed = [35:10:65];
r = [1000:200:1600];

tic
for i = 1:numel(speed)
    for j = 1:numel(r)
        calc(i,j) = batterysize(cell_2170,empty_mass,speed(i),r(j),crr,p);
        num_cells(i,j) = calc(i,j).num_cells(end);
    end
end

figure(1) 
h = carpet(speed, r, num_cells, 3);
Ytick = strsplit(num2str(get(gca, 'Ytick')));
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
ylabel('Number of Cells')
box on

% dry_mass = [1200:200:2400]*0.453592;
% 
% for j = 1:numel(dry_mass)
%     calc(j) = batterysize(cell_2170,dry_mass(j),v,range,crr,p);
%     num_cells(j) = calc(j).num_cells(end);
% end
% 
% 
% figure(2) 
% plot(dry_mass/0.453592,num_cells)
% xlabel('Empty Weight (lbs)')
% ylabel('Number of Cells')
% box on

num_cells = [4000:1000:10000];
speed = linspace(35,65,7);

for i = 1:numel(speed)
    for j = 1:numel(num_cells)
        range(i,j) = rangeCalc(num_cells(j),cell_2170,empty_mass,speed(i),crr,p);
    end
end

figure(3) 
h = carpet(num_cells,speed, range, 3);
Ytick = strsplit(num2str(get(gca, 'Ytick')));
curtick = get(gca, 'YTick');
set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));
ylabel('Range (miles)')
box on
toc

%% Functions

function calc = batterysize(cell,empty_mass,v,range,crr,p)
g = 9.81;
rho = 1.225; 
v = v*.44704;
t = range/(v/.44704); 
num_stops = round(t/4);

cfd_v = 65*.44704; 
drag = 190; 

CdA = (2*drag)/(rho*cfd_v^2);

i = 1; 

% First iteration calculates empty car values                                                       
calc.m(i) = empty_mass;                              % Sets mass to empty car weight [kg]

calc.num_cells(i) = 0;

calc.Fd(i) = .5*rho*CdA*v^2 + crr*g*calc.m(i)*cosd(2.9);                 % force required to overcome drag [N]

calc.P(i) = (calc.Fd(i)*v)/1000;                   % power (kW)
calc.P(i) = calc.P(i)/.85;                           % factor in powertrain losses

calc.kwh(i) = calc.P(i)*t;                           % capacity
calc.need(i) = calc.kwh(i);                          % how much extra capacity needed

while calc.need(i) > eps
    i = i+1;
    % add a cell to subtract energy needed
    calc.m(i) = calc.m(1) + cell.mass/1000*(i-1); % add weight of cell to total mass (kg)
    
    calc.num_cells(i) = calc.num_cells(i-1)+1; % add to cell counter
    
    calc.Fd(i) = .5*rho*CdA*v^2 + crr*g*(calc.m(i)+DFcalc(v)/g)*cosd(2.9); % force required to overcome drag (Newtons)
    calc.P(i) = (calc.Fd(i)*v)/1000; % power (kW)
    calc.P(i) = calc.P(i)/.85; % factor in powertrain losses
    acc_kwh(i) = AccelerationPower(v,calc.m(i),CdA,p);
    calc.kwh(i) = calc.P(i)*t + num_stops*acc_kwh(i); % capacity
    calc.need(i) = calc.kwh(i) - cell.energy/1000*(i-1); % how much extra capacity needed    
    
end


end

%%
function T_out = TfromRPM(input_RPM,p)

RPM = [64 110 134 146 160 171 183 193 202 211 221 233 249 275 285 295 308 323 335 ...
    350, 366 386, 405, 425, 445, 468, 494, 516, 545, 570, 592, 627, 672, 711, ...
    766, 835, 856, 862, 870, 875, 879, 883, 887, 892, 896, 900, 902, 904, 909, ...
    913, 918, 923, 925, 928, 930, 932, 933];


T = [279.6, 278, 275.5, 273.2, 270.8, 268, 265.8, 264.2, 261.5, 258.7, 255.7, 252.6 ...
    248.2, 239.8, 235.7, 231.6, 227.2, 221.8, 216.9, 212.1, 207.3, 200.8, 195.1 ...
    189.5, 183.6, 178, 171.3, 165.6, 158.5, 152.2, 147.8, 142.9, 136, 129.7, 122.8 ...
    114.7, 107.5, 100.7, 91.2, 83.2, 77, 70.1, 64.2, 56.9, 51.4, 46, 40, 34.3, 29.2 ...
    23.9, 18.5, 11.8, 7.9, 5.5, 3.4, 2.3, .7];


% Find y-intercept
p2 = polyfit(RPM(1:2),T(1:2),1);
y = polyval(p2,0);

RPM = [0 RPM];
T = [y T];

i = 1;
ref_RPM = RPM(1);

while input_RPM > ref_RPM
    i = i+1;
    ref_RPM = RPM(i);
end

if input_RPM == ref_RPM
    T_out = polyval(p(i,1:end),input_RPM);
else
    T_out = polyval(p(i-1,1:end),input_RPM);
end

end
%%
function [kwh] = AccelerationPower(vmax,m,CdA,p)

rho = 1.225; 
tireDiameter = 27.6*0.0254; % [m]
v = 0;
x = 0;
RPM = 0;
g = 9.81;
t_step = 1; % s
P = 0;
i = 1; 
Fnet = 0;
drag = 0;
crr = .0068;
a = 0;
t = 0;
T = 0;

vRPM = vmax*60/(pi*tireDiameter);

while RPM(i) < vRPM

    T(i) = 2*TfromRPM(RPM(i),p);
    
    Fwheel(i) = T(i)/(tireDiameter/2); % force at wheel
    drag(i) = .5*rho*CdA*v(i)^2;
    rollResis = crr*(m+DFcalc(v(i))/g)*g; 
    Fnet(i) = Fwheel(i) - drag(i) - rollResis;

    P(i+1) = Fnet(i)*v(i);
    
    t(i+1) = t(i) + t_step; 
    a(i+1) = Fnet(i)/m;
    v(i+1) = v(i) + a(i)*t_step;

    RPM(i+1) = v(i)*60/(pi*tireDiameter);
    i = i+1;
end

P = P(end)/.85;
time = t(end)/3600;

kwh = P*time/1000;
end
%%
function range = rangeCalc(num_cells,cell,empty_weight,v,crr,p)
    v = v*.44704; 
    rho = 1.225; 
    g = 9.81; 
    total_mass = num_cells*cell.mass/1000 + empty_weight; 
    
    cfd_v = 65*.44704; 
    cfd_drag = 190; 

    CdA = (2*cfd_drag)/(rho*cfd_v^2);
   
    
    kwh = num_cells*cell.energy/1000; 
    drag = .5*rho*CdA*v^2;
    roll_res = (total_mass+DFcalc(v)/g)*g*crr;
    
    F = drag + roll_res;
    
    power = (F*v)/1000; % power (kW)
    power = power/.85; % factor in powertrain losses
    
    t = kwh/power;
    
    num_stops = round(t/4); 
    
    kwh = kwh - num_stops*AccelerationPower(v,total_mass,CdA,p);
    
    t = kwh/power*3600; 
    
    range = v*t*0.000621371;
end
    
%% 
function downforce = DFcalc(v) 
    DF = [0 1100];
    vel = [0 65*.44704];
    
    p = polyfit(vel,DF,2);
    downforce = polyval(p,v);
end


