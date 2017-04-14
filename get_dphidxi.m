function [ dphidxi ] = get_dphidxi( gauss_points, num_phis )
%get_dphidxi Returns values of dphidxi at gauss points
%   Returns values in an array where row number is the phi, and column
%   number is the gauss point

% Create a table of derivatives of phi
% Rows: phi1 through phi8
% columns gauss points 1 - 8
% take the value gauss_points(i,2) for the eta component
num_gp = size(gauss_points,1);
dphidxi = zeros(num_phis,num_gp);
xicount = 1;
for i =1:2
    for j=1:2
        for k=1:2
            for l=1:num_gp
                dphidxi(xicount,l) = (-1)^(i)*0.5*(1+(-1)^(k)...
                *gauss_points(l,2))*0.5*(1+(-1)^(j)*gauss_points(l,3))*0.5; 
                
            end
            xicount = xicount + 1;
        end
    end
end

end

