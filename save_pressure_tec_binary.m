% This file will convert csv data files from my simulations into tecplot
% format
% xdims number of nodes in x-direction
% ydims number of nodes in y-direction
% flag is to test gmsh import, if true then no other variables
% flux is true when data being passed are IgA concentrations with IgA flux
% also calculated and passed to the function

function save_pressure_tec_binary(xpoints, ypoints, zpoints,...
    xdims, ydims, zdims, pressure, flag, flux, varargin)
cd('/home/sarah/Downloads/Tecplot conversion');
% Folder where data file will be saved
savefolder = '/home/sarah/Documents/Research Writing/Counter_current_paper/';
%cd(savefolder);
nx = xdims;
ny = ydims;
nz = zdims;
Nnodes = length(xpoints); % total number of nodes

Nels = (nx-1)*(ny-1)*(nz-1); % Number of elements in domain
tdata = [];
extraArgs = length(varargin);
if extraArgs ~= 1 && ~flag 
    Nvar = 4 + extraArgs;
elseif flag
    Nvar  = 4;
else
    Nvar = 4; % Number of variables - pressure and x,y,z info
end

tdata.Nvar = Nvar;
if Nvar == 4
    tdata.varNames = {'X', 'Y', 'Z', 'IgA'}; % Load the variable names
elseif extraArgs == 3 && ~flux
    tdata.varNames = {'X', 'Y', 'Z', 'Pressure', 'U', 'V', 'W'};
elseif extraArgs == 3 && flux
    tdata.varNames = {'X', 'Y', 'Z', 'IgA', 'Fx', 'Fy', 'Fz'};
elseif flag
    tdata.varNames = {'X','Y','Z','Test'};
else
    error('Error. Unrecognized number of input variables');
end

tdata.vformat(1:Nvar) = 1; %Sets variables to float

if ~flag
    tdata.FEvolumes.varlock = 0; % Sets variable to be nodal (1 gives cell centers)
    tdata.FEvolumes.x = xpoints;
    tdata.FEvolumes.y = ypoints;
    tdata.FEvolumes.z = zpoints;
else
    tdata.FEsurfaces(1).varlock = 0; % Nodal data
    tdata.FEsurfaces(1).x = xpoints;
    tdata.FEsurfaces(1).y = ypoints;
    tdata.FEsurfaces(1).z = zpoints;
end

if extraArgs ~= 0 && ~flag
    variableData = [pressure  varargin{1} varargin{2} varargin{3}];
    tdata.FEvolumes.v = variableData';
elseif flag
    tdata.FEsurfaces(1).e2n = varargin{1};
    tdata.FEsurfaces(1).v = ones(1, Nnodes);
    tdata.FEsurfaces(1).order = 4;
else
    tdata.FEvolumes.v = reshape(pressure, [1, Nnodes]); % Load the pressure data as a variable
end

if ~flag
    % Create the connectivity list that defines elements
    connect_list = zeros(Nels,8);
    element = 1;
    for elz = 1:zdims-1
        for ely=1:ydims-1
            for elx=1:xdims-1
                connect_list(element,:) = [ydims*(elx-1)+1+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*(elx-1)+2+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*nx+ydims*(elx-1)+2+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*nx+ydims*(elx-1)+1+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*elx+1+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*elx+2+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*nx+ydims*elx+2+(ely-1)+(elz-1)*ydims*nx, ...
                    ydims*nx+ydims*elx+1+(ely-1)+(elz-1)*ydims*nx];
                element = element +1;
            end
        end
    end

    
    tdata.FEvolumes.e2n = connect_list;
end

filename = '/home/sarah/Documents/Research Writing/Counter_current_paper/bm04_phi_05_r_98_conv_flux.plt';

mat2tecplot(tdata, filename);
end