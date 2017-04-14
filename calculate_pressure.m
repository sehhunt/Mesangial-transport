function [P] = calculate_pressure( Pmatrix, numx_nodes,...
    numy_nodes,numz_nodes, ...
    node_coord_matrix, kappa, mu)
% Future work - add a num_gp parameter that determines which Gauss points
% to use
%calculate_pressure This function calculates the stiffness matrix for the
%Darcy equation
%   Detailed explanation goes here

NOP_array = NOP(numy_nodes, numz_nodes,numx_nodes);
NNODES = numy_nodes*numx_nodes*numz_nodes;
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


 element = 1;
for elx=1:numx_nodes-1
    for ely = 1:numy_nodes-1
        for elz =1:numz_nodes-1
            
            % Now we loop over the gauss points. We are evaluating the integral
            % using gauss quadrature, so we sum contributions to the integral at
            % each gauss point
            for gp = 1:num_gp
                % Get the xi, eta, zeta coordinates of this gauss point
                %gp_coords = gauss_points(gp,1:2);
                gp_weight = gauss_points(gp, 4); % column 4 of gauss_points is gp weight
                dphidxis = dphidxi(:,gp); % get every value of dphidxi at this gp
                % This is a column vector
                dphisdetas = dphideta(:,gp); %similar to above
                dphisdzetas = dphidzeta(:,gp);
                
                el_nodes = NOP_array(element,:); 
                % This returns a length 8 vector with the global node numbers
                % of the element elx, ely, elz. The first element of the vector
                %is the global node number of the first node for the element, etc.
            
                % Column 1 of node_coords is x coordinate of node # row #
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
                % Calculate dphi/dx dphi/dy etc using Jacobian, for each
                % phi (will end up with eight) These are stored in a matrix
                % where each row is the phi, e.g., row 1 is derivatives of
                % phi 1, row 2 is derivatives of phi 2 etc. Columns are
                % derivatives with respect to x,y,z in column 1,2,3
                dphi = zeros(8,3);
                for ph = 1:num_gp
                    dphi(ph, :) = Jacobian2\[dphidxis(ph); dphisdetas(ph); dphisdzetas(ph)];
                end
                                
                
                % Get all the phis that we iterate over 
                for i=1:num_phis % here num_phis substitutes for number of
                    % local nodes
                    for j=1:num_phis
                        Pmatrix(el_nodes(i), el_nodes(j)) = ...
                            Pmatrix(el_nodes(i), el_nodes(j))  -...
                            kappa(element)/mu*gp_weight*abs(J)*(dphi(i,1)*dphi(j,1)+...
                            dphi(i,2)*dphi(j,2)+dphi(i,3)*...
                            dphi(j,3));
                    end
                end
                
            end
            element = element + 1;
        end
   end
end

P = Pmatrix;
end

