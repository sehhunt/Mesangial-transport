function [ dQdz ] = calculate_capillary_Q_no_mes(k,R,Lz,Q0,m,direction,c0,...
    P_cap,P_urinary,P_max)
%calculate_capillary_Q_no_mes Returns dQdz in the case of no mesangium.
%Output is non-dimensional
%
%   k: permeability of capillary wall
%   R: radius of capillary
%   Lz: axial length of capillary
%   Q0: volumetric flow in capillary at entrance
%   m: width of mesangial area abutting capillary
%   direction: +1 means flow in + z direction, -1 means flow in negative
%   z-direction
%   P_cap: pressure in capillary as a function of z, in Pascals
%   P_urinary: pressure in urinary space, in mmHg
%   P_max: pressure that is set to zero, in mmHg

a1 = 1.63;
a2 = 0.294;

dQdz = @(conc, z_pos) (-direction*2*k*R*Lz*P_max*133.3/Q0*...
    (pi-asin(m/(2*R)))*((P_cap(z_pos)/133.3-P_urinary)/P_max-...
    a1*c0*conc/P_max-a2*c0^2*conc^2/P_max));

end

