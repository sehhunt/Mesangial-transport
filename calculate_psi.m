function [ psi ] = calculate_psi( v_bm, v_m, delta_bm, mes_length)
%calculate_psi Calculate the dimensionless number that compares the
%resistance of the mesangial matrix to the basement membrane
%   psi = kappa_bm * mes_length/(kappa_m * delta_bm)
%   v_bm = volume fraction of basement membrane
%   v_m = volume fraction of mesangial matrix
%   delta_bm - thickness of basement membrane
%   mes_length - x length of mesangial matrix
%   mu - viscosity of blood plasma.

r_fibers = 2*10^(-9); % fiber radius, meters
r_fibers_bm = 2*10^(-10); % BM fiber radius, meters

km = 3*r_fibers^2/(20*v_m)*(-log(v_m)...
    -0.931)*10^(12); % um^2
kbm = 3*r_fibers_bm^2/(20*v_bm)*(-log(v_bm)...
    -0.931)*10^(12); % um^2

psi = kbm*mes_length/(km*delta_bm);
end

