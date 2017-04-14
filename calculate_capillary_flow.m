function [Q ] = calculate_capillary_flow( dQdz, capillary_conc,dz, nnodes,...
    direction)
%calculate_capillary_flow Summary of this function goes here
%   Detailed explanation goes here
Q = zeros(nnodes,1);
if direction > 0
    Q(1) = 1; % dimensionless value
    for i =2:nnodes
        Q(i) = Q(i-1)+...
            dz*dQdz(capillary_conc(i-1),i-1);
    end
else
    Q(nnodes) = 1;
    for i=nnodes-1:-1:1
        Q(i) = Q(i+1)+...
            dz*dQdz(capillary_conc(i+1),i+1);
    end
end

end

