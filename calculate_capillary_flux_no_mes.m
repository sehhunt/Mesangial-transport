function [ dQ ] = calculate_capillary_flux_no_mes( direction, k, R, m,...
    Lz, P_cap, P_urinary, P_max, Q0)
%calculate_capillary_flux_no_mes Calculate the flux across the capillary
%wall with no mesangium present. Returns a function dQdz of concentration
%and pressure
%   direction - whether the flow is in the positive z-direction (+1) or the
%   negative z-direction (-1)
%   k permeability of the capillary wall
%   R radius of the capillary
%   Lz length of the capillary
%   P_max reference pressure
%   P_urinary Pressure in the urinary space, normalized to P_max
%   Q0 entrance volumetri flow
%   
%   This code assumes that P_urinary and P_max are in mmHg, and that the
%   P_cap vector is in Pascals. This results in the random factors of
%   133.3, which convert the Pascals, to mmHg. It needs to calculate in
%   mmHg, because the formula for osmotic pressure, (a1*c0*conc +
%   a2*c0^2*conc^2), outputs pressure in mmHg

dQ = @(conc, z_pos) (-direction*2*k*R*Lz*P_max*133.3/Q0*(pi-asin(m/(2*R)))...
    *((P_cap(z_pos)*1/133.3-P_urinary)/P_max-...
    ((a1*c0*conc+a2*c0^2*conc^2))/P_max));
end

