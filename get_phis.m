function [phis] = get_phis(gauss_points,num_phis)
% phis Returns the value of the phis function (linear basis functions) at
% the given gauss points.
% Row: which phi we have, e.i. phi1, phi2, etc.
% Column: what gausspoint we have
% Create a table of values of phi
num_gp = size(gauss_points,1);
phis = zeros(num_phis,num_gp);
phicount = 1;
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:num_gp
                phis(phicount,l) = (1+(-1)^(i)*gauss_points(l,1))/2*...
                (1+(-1)^(k)*gauss_points(l,2))/2*(1+(-1)^(j)*...
                gauss_points(l,3))/2;
            end
            phicount = phicount + 1;
        end
    end
end

end