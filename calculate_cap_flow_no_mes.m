function [ capillary_flow] = calculate_cap_flow_no_mes(...
    dz, flow_func, nnodes, conc_func)
%calculate_cap_flow_no_mes This function uses a forward Euler
%method to calculate the capillary concentration for an equation of the
%form dc/dz = f(c_cap, z_pos)
%   inlet_concentration is the incoming concentration of the protein of
%   interest
%   
%   dz is the z spacing of the nodes. This should be a dimensionless number
%
%   cfunc is the f(c_cap, z_pos) mentioned above. This is the function that
%   determines the total capillary concentration. It needs to be a function
%   of two variables, the first is the concentration, and the second is the
%   z position, numbered from 1 at the z-low end to nnodes at the z-high
%   end. This function should return the dimensionless concentration
%
%   nnodes is the number of nodes in the z-direction that we need to fill
%   in
capillary_flow = zeros(length(nnodes),1);
capillary_flow(1) = 1; % no dimensions!
for i =2:nnodes
    capillary_flow(i) = capillary_flow(i-1)+...
        dz*flow_func(conc_func(i-1),i-1);
end

end