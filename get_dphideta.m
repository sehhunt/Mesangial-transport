function [ dphideta ] = get_dphideta( gauss_points, num_phis)
%get_dphideta Returns the dphidetas 
%   Detailed explanation goes here
% Make a table of dphi/deta
num_gp = size(gauss_points,1);
dphideta = zeros(num_phis, num_gp);
etacount = 1;
for i=1:2
    for j=1:2
        for k = 1:2
            for l=1:num_gp
                dphideta(etacount,l) = 0.5*(-1)^(k)*...
                    (1+(-1)^(i)*gauss_points(l,1))*0.5*(1+(-1)^(j)*...
                    gauss_points(l,3))*0.5;
                
            end
            etacount= etacount+1;
        end
    end
end

end

