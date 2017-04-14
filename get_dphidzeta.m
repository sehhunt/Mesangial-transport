function [ dphidzeta ] = get_dphidzeta( gauss_points, num_phis )
%get_dphideta Returns values of dphideta at gauss points
%   Returns values in an array where row number is the phi, and column
%   number is the gauss point

% Create a table of derivatives of phi
% Rows: phi1 through phi8
% columns gauss points 1 - 8
% take the value gauss_points(i,2) for the eta component

% Create a table of dphi/dzeta
num_gp = size(gauss_points,1);
dphidzeta = zeros(num_phis, num_gp);
zetacount = 1;
for i=1:2
    for j=1:2
        for k=1:2
            for l =1:num_gp
                dphidzeta(zetacount, l) = 0.5*(-1)^(j)*(1+(-1)^(i)*...
                    gauss_points(l,1))*0.5*(1+(-1)^(k)*gauss_points(l,2))*...
                    0.5;
            end
            zetacount = zetacount + 1;
        end
    end
end

end
