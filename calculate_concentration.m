function [C] = calculate_concentration(Pvector,Cmatrix, numx_nodes,...
    numy_nodes,numz_nodes,...
     node_coord_matrix,kappa,mu, D) %#codegen
% Calculates the concentration distrubution, given the pressure
% distribution
%   calculate_concentration is the loop that calculates the stiffness
%   matrix? (the A matrix) for the convection diuffusion equation

% Get the array that has global node numbers:
NOP_array = NOP(numy_nodes,numz_nodes,numx_nodes);
num_gp = 8;
gauss_points = zeros(8,4);
    gauss_points(1,:) = [-0.5774 -0.5774, -0.5774  1];
    gauss_points(2,:) = [-0.5774 0.5774 -0.5774 1];
    gauss_points(3,:) = [-0.5774 -0.5774 0.5774 1];
    gauss_points(4,:) = [-0.5774 0.5774  0.5774 1];
    gauss_points(5,:) = [0.5774 -0.5774, -0.5774  1];
    gauss_points(6,:) = [0.5774 0.5774, -0.5774  1];
    gauss_points(7,:) = [0.5774 -0.5774, 0.5774  1];
    gauss_points(8,:) = [0.5774 0.5774, 0.5774  1];
  
num_phis = 8; % Also the number of local nodes

% Create a table of derivatives of phi
% Rows: phi1 through phi8
% columns gauss points 1 - 8
% take the value gauss_points(i,2) for the eta component
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

% Make a table of dphi/deta
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
% Create a table of dphi/dzeta
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
% Create a table of values of phi
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

% Loop 1: loop over elements (this actually takes three loops as we loop in
% x, y, and z
element = 1;
for elx = 1:numx_nodes-1
    for ely=1:numy_nodes-1
        for elz=1:numz_nodes-1
                % This returns a length 8 vector with the global node numbers
                % of the element elx, ely, elz. The first element of the vector
                %is the global node number of the first node for the element, etc.
            el_nodes = NOP_array(element,:);
            
               for gp =1:num_gp
                % The first part of this loop is the same because we still
                % have a dphi_i dotted to a dphi_j
                gp_weight = gauss_points(gp, 4); % column 4 of gauss_points is gp weight
                dphidxis = dphidxi(:,gp); % get every value of dphidxi at this gp
                % This is a column vector
                dphisdetas = dphideta(:,gp); %similar to above
                dphisdzetas = dphidzeta(:,gp);
                
               
                % Calculate components of Jacobian
                
                dxdxi= sum(node_coord_matrix(el_nodes,1).*dphidxis);
                dxdeta = sum(node_coord_matrix(el_nodes,1).*dphisdetas);
                dxdzeta = sum(node_coord_matrix(el_nodes,1).*dphisdzetas);
                dydxi = sum(node_coord_matrix(el_nodes,2).*dphidxis);
                dydeta = sum(node_coord_matrix(el_nodes,2).*dphisdetas);
                dydzeta = sum(node_coord_matrix(el_nodes,2).*dphisdzetas);
                dzdxi = sum(node_coord_matrix(el_nodes,3).*dphidxis);
                dzdeta = sum(node_coord_matrix(el_nodes,3).*dphisdetas);
                dzdzeta = sum(node_coord_matrix(el_nodes,3).*dphisdzetas);
                Jacobian = [dxdxi dxdeta dxdzeta; dydxi dydeta dydzeta;...
                    dzdxi dzdeta dzdzeta];
                Jacobian2 = [dxdxi dydxi dzdxi; dxdeta dydeta dzdeta;...
                    dxdzeta dydzeta dzdzeta];
                J = det(Jacobian);
                
                phi_at_gp = phis(:,gp); % This returns all the phis (1 - 8)
                % at this particular gauss point
                
                              
                
                % Calculate dphi/dx dphi/dy etc using Jacobian, for each
                % phi (will end up with eight) These are stored in a matrix
                % where each row is the phi, e.g., row 1 is derivatives of
                % phi 1, row 2 is derivatives of phi 2 etc. Columns are
                % derivatives with respect to x,y,z in column 1,2,3
                dphi = zeros(8,3);
                for ph = 1:num_gp
                    dphi(ph, :) = Jacobian2\[dphidxis(ph); dphisdetas(ph); dphisdzetas(ph)];
                end
                % Now, sum the pressure nodes with the dphis to get the
                % gradient of the pressure at the gauss point gp
                
                pressures = Pvector(el_nodes);
                grad_p_x = sum(dphi(:,1).*pressures);
                grad_p_y = sum(dphi(:,2).*pressures);
                grad_p_z = sum(dphi(:,3).*pressures);
                % Now construct the matrix
                for i=1:num_phis % here num_phis substitutes for number of
                    % local nodes
                   for j=1:num_phis
                        Cmatrix(el_nodes(i), el_nodes(j)) = ...
                            Cmatrix(el_nodes(i), el_nodes(j)) -...
                            gp_weight*abs(J)*D*(dphi(i,1)*dphi(j,1)+... % should this be abs(J)?
                            dphi(i,2)*dphi(j,2)+dphi(i,3)*...
                            dphi(j,3)) -...
                            gp_weight*abs(J)*kappa(element)/mu*(dphi(i,1)*phi_at_gp(j)*grad_p_x +...
                            dphi(i,2)*phi_at_gp(j)*grad_p_y+dphi(i,3)*phi_at_gp(j)...
                            *grad_p_z);
                    end
                end
               end
            
            element = element + 1;
        end
    end
end

C = Cmatrix;
end

