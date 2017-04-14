function [ dQ ] = with_mes_concentration( R, k, h, Lz, Q0, c0, m, P_cap, ...
    plasma_func, P_urinary, P_max, direction)
%with_mes_concentration Calculates dQ/dz as a function of z-position based
%on the loss of plasma across the capillary wall as well as flux from the
%mesangium. It is formatted to include the mesangium, that is it returns 
%dQ/dz, not dc/dz
%   plasma_func is a function of z-position that gives total flux into or
%   out of capillary
% direction tells the direction of flow for this capillary. +1 is flow in
% the positive z-direction, and -1 is flow in the negative z-direction
% 

a1 = 1.63;
a2 = 0.294;


dQ = @(conc, z_pos) (direction*2*h*Lz/Q0*plasma_func(z_pos)-...
    direction*2*k*R*Lz*P_max*133.3/Q0*(pi-asin(m/(2*R)))...
    *((P_cap(z_pos)*1/133.3-P_urinary)/P_max-...
    ((a1*c0*conc+a2*c0^2*conc^2))/P_max));
% extra factor of 1/133.3 to convert P_cap, in Pascals, to mmHg

end

