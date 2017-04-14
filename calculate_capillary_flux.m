function [ dQcdz ] = calculate_capillary_flux( Q0, c0, Lz, h, ...
    conc_flux )
%calculate_capillary_flux This function creates a function for dQc/dz From
%the mesangial surface concentration flux
%   Detailed explanation goes here

dQcdz = @(z_pos) (Lz/(c0*Q0)*2*h*conc_flux(z_pos));
end

