function [ capillary_concentration ] = calculate_cap_concentration_no_mes(...
    dz, cfunc, nnodes, direction, varargin)
%calculate_cap_concentration_no_mes This function uses a forward Euler
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
%   direction tells whether flow is in positive (+1) or negative (-1)
%   z-direction
%
%   nnodes is the number of nodes in the z-direction that we need to fill
%   in
if length(varargin) ==1
    b = varargin{1}; % This is the concentration of albumin, DIMENSIONAL, in g/dL
    capillary_concentration = zeros(length(nnodes),1);
    if direction > 0
        capillary_concentration(1) = 1; % dimensionless value
        for i =2:nnodes
            capillary_concentration(i) = capillary_concentration(i-1)+...
                dz*cfunc(capillary_concentration(i-1),b(i-1),i-1);
        end
    else
        capillary_concentration(nnodes) = 1;
        for i=nnodes-1:-1:1
            capillary_concentration(i) = capillary_concentration(i+1)-...
                dz*cfunc(capillary_concentration(i+1),b(i+1),i+1);
        end
    end
else
    capillary_concentration = zeros(length(nnodes),1);
    if direction > 0
        capillary_concentration(1) = 1; % dimensionless value
        for i =2:nnodes
            capillary_concentration(i) = capillary_concentration(i-1)+...
                dz*cfunc(capillary_concentration(i-1),i-1);
        end
    else
        capillary_concentration(nnodes) = 1;
        for i=nnodes-1:-1:1
            capillary_concentration(i) = capillary_concentration(i+1)-...
                dz*cfunc(capillary_concentration(i+1),i+1);
        end
    end
end



end

