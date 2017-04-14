function [ concentration ] = no_mes_concentration( R, k, Q0, c0, m, Lz, P_cap, ...
    P_urinary, P_max, direction, varargin)
%no_mes_concentration Return the function equal to dc/dz for a capillary
%with no mesangial flow
%   R, radius of capillary.   k permeability of capillary wall, Q0 volumetric
%   flow rate at the entrance (this is constant without a
%   mesangium), m is the width of the total mesangial channel P_ cap
%   capillary pressure (function of z). P_cap is a matrix, giving the
%   capillary pressure at the node in position i. That is P_cap(i) is the
%   pressure in the capillary at node i (counting in z, from low to high).
%   It returns a function, concentration, that gives dc/dz as a function of
%   capillary concentration (conc) and position along the length of the
%   capillary (z_pos). conc must be in g/100 mL. P_cap must be in Pascals
%   P_urinary is the tubule pressure. It must be in mmHg. Lz is the length
%   in the z-direction, to non-dimensionalize z. c0 is the inlet flow
%   concentration. P_max is in mmHg

% All outputs are non-dimensional

a1 = 1.63;
a2 = 0.294;

if length(varargin) == 1
   % This is for IgA, or any protein that doesn't contribute to the osmotic
   % pressure
    concentration = @(conc, b, z_pos) (direction*2*k*R*Lz*P_max*133.3/Q0*(pi-asin(m/(2*R)))*conc^2 ...
        *((P_cap(z_pos)*1/133.3-P_urinary)/P_max-((a1*b+...
        a2*b^2))/P_max));
    % extra factor of 1/133.3 to convert P_cap, in Pascals, to mmHg
    
else
    concentration = @(conc, z_pos) (direction*2*k*R*Lz*P_max*133.3/Q0*(pi-asin(m/(2*R)))*conc^2 ...
        *((P_cap(z_pos)*1/133.3-P_urinary)/P_max-((a1*c0*conc+...
        a2*c0^2*conc^2))/P_max));
    % extra factor of 1/133.3 to convert P_cap, in Pascals, to mmHg
    
end

end
