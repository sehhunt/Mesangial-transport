% This is the main function for running a simulation of transport within
% (for now) the mesangial slab

% Clear variables
%clear all;

 % Change the matlab path to the right directory to run convection
 % diffusion code
 
 %cd('/panfs/roc/groups/13/barocasv/hunts/Documents/mesangial_flow/trunk');
 cd('/home/barocasv/hunt0646/Documents/mesangial_flow/trunk')
addpath('VTK conversions/')
% The first thing to do is to create the mesh

% How thick of a basement membrane do you want?

bm_thickness = 0.4; % in um

% How wide is the mesangial channel?
mes_thickness = 0.5;

% The basement membrane is the top in the y direction, the code below makes
% it begin at y = 0.5 
% No longer have explicit basement membrane
yvalues = [0:0.05:mes_thickness];
ylimits = [0 mes_thickness];

% the x direction is the length of the mesangium towards the interior of
% the glomerulus

xvalues = [0:0.1:5]; % um 
xlimits = [0 5];
% the z direction is along the glomerular capillary
Lz = 20;

zvalues  = 0:0.5:Lz; % um
zlimits = [0 Lz];

capillary_zs = 0:0.005:Lz;
cap_dz = 0.005;
num_cap_zs = length(capillary_zs);

% Create a mesh of squares
[x,y,z] = meshgrid(xvalues,yvalues,zvalues);

% Get number of nodes in each direction
x_nodes = size(x,2);
y_nodes = size(y,1);
z_nodes = size(z,3);

% Get total number of nodes
NUM_NODES = x_nodes*y_nodes*z_nodes;
NUM_ELEMENTS = (x_nodes-1)*(y_nodes-1)*(z_nodes-1);
% Initialize the permeability and viscosity parameters for the given mesh
% Initialize the permeability and viscosity parameters for the given mesh
mu = 1.35 * 10^(-3); %Pa s 1.25*(10^-9); % kg/(s um) apparant viscosity of blood plasma

% Calculate D and permeability using fiber radius, volume fraction, and
% solute size
rIgA = 98.25*10^(-10); % This needs to be in meters
ralbumin = 13.7*10^(-10); % This needs to be in meters
T = 310; % Kelvin
kb = 1.3806*10^(-23); %m^2 kg/(s^2 K)
D0_alb = kb*T/(6*pi*mu*ralbumin)*10^12; % um^2/s
D0_IgA = kb*T/(6*pi*mu*rIgA)*10^12; % um^2/s


volume_fraction_bm = 0.3;
volume_fraction_matrix = 0.05;
r_fibers = 2*10^(-9); % fiber radius, meters
r_fibers_bm = 2*10^(-10); % BM fiber radius, meters
km = 3*r_fibers^2/(20*volume_fraction_matrix)*(-log(volume_fraction_matrix)...
    -0.931)*10^(12); % um^2
kbm = 3*r_fibers_bm^2/(20*volume_fraction_bm)*(-log(volume_fraction_bm)...
    -0.931)*10^(12); % um^2
lambda_IgA = rIgA/r_fibers; % unitless
lambda_albumin = ralbumin/r_fibers; %unitless

b_albumin = 0.174*log(59.6/lambda_albumin);
b_IgA = 0.174*log(59.6/lambda_IgA);
f_albumin = (1+lambda_albumin)^2*volume_fraction_matrix;
f_IgA = (1+lambda_IgA)^2*volume_fraction_matrix;
D_albumin = D0_alb*exp(-pi*volume_fraction_matrix^b_albumin)*...
    exp(-0.84*f_albumin^1.09);
D_IgA = D0_IgA*exp(-pi*volume_fraction_matrix^b_IgA)*...
    exp(-0.84*f_IgA^1.09);
% kappa(:) = ones(NUM_ELEMENTS,1).*km;
% 
% element_nodes = NOP(y_nodes,x_nodes,z_nodes);
% 
% for elx=1:x_nodes-1
%     for ely=1:y_nodes-1
%         for elz=1:z_nodes-1
%             coords=element_nodes(elx,ely,elz,num_nodesy,num_nodesz);
%             if node_coords(coords(2),2)>0.5 % This assumes that the basement
%                 % membrane begins at y = 0.5 --> This should be a PARAMETER
%                 % CHANGE
%                 kappa(element) = kbm;
%             end
%             element = element+1;
%         end
%     end
% end

% Initialize the Pressure and concentration values
% P = zeros(NUM_NODES, 1);
% C = zeros(NUM_NODES, 1);
%%
direction_left = +1;
direction_right = -1;

reference_pressure = 35.3; %mmHg
epsilon = 0.05;
pressure_function = @(z)(reference_pressure+direction_left*epsilon*0.5*reference_pressure-...
    direction_left*epsilon*reference_pressure*z/Lz);

reference_pressure2 = reference_pressure;
right_pressure = @(z)(reference_pressure2+direction_right*epsilon*1.5*reference_pressure2-...
    direction_right*epsilon*reference_pressure2*z/Lz);

right_radius = 5; %um
left_radius = 5; %um
m=3; % um

area_left = 2*(pi-asin(m/(2*left_radius)))*left_radius*Lz; % um^2
area_right = 2*(pi-asin(m/(2*right_radius)))*right_radius*Lz;
h = mes_thickness; % um
Q0L = 8.3*10^5; %um^3/s
Q0R = 8.3*10^5; %um^3/s



k = 5*10^(-2);%F*Q0L/(area_left*reference_pressure*133.3); %um/(Pa s);
kR = k;
F = k*area_left*reference_pressure*133.3/Q0L;

flag = true;
cR = 1;%5.8; %g/100 mL
cL = 5.8;

% Need to add an IgA concentration to these simulations
cIgA_R = 1; % g/100 mL
cIgA_L = 0.2;

issingle = false;
%% Run convection-diffusion code
% Now that we've created the mesh and specified the boundary conditions, we
% can solve for the pressure and concentration


[QL, QR, P, C, CIg, alb_capL, alb_capR,loop_counter, outer_count] = cross_capillary(xvalues,...
    yvalues,zvalues,pressure_function,right_pressure,left_radius,right_radius,...
    k,kR,m,h, Q0R,cR,Q0L,cL,cIgA_R, cIgA_L,direction_left,direction_right, ...
    bm_thickness,mu,km,kbm,D_albumin,D_IgA,capillary_zs, cap_dz, issingle, flag);
%% Output the results to a file

save_PressureVTK_binary(P,x,y,z,'/panfs/roc/groups/13/barocasv/hunt0646/Documents/Counter_sims/04012017/bm_04_r_98_phi_05_pressure.vtk');
save_ConcentrationVTK_binary(C,x,y,z,'albumin','/panfs/roc/groups/13/barocasv/hunt0646/Documents/Counter_sims/04012017/bm_04_r_98_phi_05_albumin.vtk');
save_ConcentrationVTK_binary(CIg,x,y,z,'IgA','/panfs/roc/groups/13/barocasv/hunt0646/Documents/Counter_sims/04012017/bm_04_r_98_phi_05_IgA.vtk');
