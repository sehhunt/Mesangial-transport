% Create a function that implements a neumann boundary condition
% boundary_node_matrix is boundaries from above - matrix that contains
% whether each node is on the boundary
function [b] = neumann_bc(boundary, boundary_coordinates, boundary_node_matrix, n_coords, value, bvalues, nny, nnx, nnz)
% Boundary is just a number telling which boundary - I need to convert that
% into a direction vector and coordinate, to use the onBoundary function

if (boundary == 1 || boundary==2)
    direction = [false true false];
elseif (boundary == 3 || boundary == 4)
    direction = [true false false];
else
    direction = [ false false true];
end

if mod(boundary,2) == 0
    boundary_coordinates = boundary_coordinates(2);
else
    boundary_coordinates = boundary_coordinates(1);
end
% In order to evaluate a Neumann boundary condition I need to create new
% gauss points on the face of the element. Rows are gauss points, columns
% are xi, eta, zeta coordinate and weight. These map to an element that goes
% from -1 to 1

% This needs to change to be a function of the boundary
gauss_points2 = zeros(4,4);
if boundary == 1 
    gauss_points2(1,:) = [-0.5774 -1 -0.5774 1];
    gauss_points2(2,:) = [-0.5774 -1 -0.5774 1];
    gauss_points2(3,:) = [-0.5774 -1 0.5574 1];
    gauss_points2(4,:) = [-0.5774 -1 0.5574 1];
    gauss_points2(5,:) = [0.5774 -1 -0.5774 1];
    gauss_points2(6,:) = [0.5774 -1 -0.5774 1];
    gauss_points2(7,:) = [0.5774 -1 0.5774 1];
    gauss_points2(8,:) = [0.5774 -1 0.5774 1];
elseif boundary == 2
   gauss_points2(1,:) = [-0.5774 1 -0.5774 1];
    gauss_points2(2,:) = [-0.5774 1 -0.5774 1];
    gauss_points2(3,:) = [-0.5774 1 0.5574 1];
    gauss_points2(4,:) = [-0.5774 1 0.5574 1];
    gauss_points2(5,:) = [0.5774 1 -0.5774 1];
    gauss_points2(6,:) = [0.5774 1 -0.5774 1];
    gauss_points2(7,:) = [0.5774 1 0.5774 1];
    gauss_points2(8,:) = [0.5774 1 0.5774 1];
elseif boundary == 3
    gauss_points2(1,:) = [-1 -0.5774 -0.5774 1];
    gauss_points2(2,:) = [-1 0.5774 -0.5574 1];
    gauss_points2(3,:) = [-1 -0.5774 0.5774 1];
    gauss_points2(4,:) = [-1 0.5774 0.5774 1];
    gauss_points2(5,:) = [-1 -0.5774 -0.5774 1];
    gauss_points2(6,:) = [-1 0.5774 -0.5574 1];
    gauss_points2(7,:) = [-1 -0.5774 0.5774 1];
    gauss_points2(8,:) = [-1 0.5774 0.5774 1];
elseif boundary == 4
    gauss_points2(1,:) = [1 -0.5774 -0.5774 1];
    gauss_points2(2,:) = [1 0.5774 -0.5574 1];
    gauss_points2(3,:) = [1 -0.5774 0.5774 1];
    gauss_points2(4,:) = [1 0.5774 0.5774 1];
    gauss_points2(5,:) = [1 -0.5774 -0.5774 1];
    gauss_points2(6,:) = [1 0.5774 -0.5574 1];
    gauss_points2(7,:) = [1 -0.5774 0.5774 1];
    gauss_points2(8,:) = [1 0.5774 0.5774 1];
elseif boundary == 5
    gauss_points2(1,:) = [-0.5774 -0.5774 -1 1];
    gauss_points2(2,:) = [-0.5774 0.5574 -1 1];
    gauss_points2(3,:) = [-0.5774 -0.5774 -1 1];
    gauss_points2(4,:) = [-0.5774 0.5574 -1 1];
    gauss_points2(5,:) = [0.5774 -0.5774 -1 1];
    gauss_points2(6,:) = [0.5774 0.5774 -1  1];
    gauss_points2(7,:) = [0.5774 -0.5774 -1 1];
    gauss_points2(8,:) = [0.5774 0.5774 -1  1];
else
    gauss_points2(1,:) = [-0.5774 -0.5774 1 1];
    gauss_points2(2,:) = [-0.5774 0.5574 1 1];
    gauss_points2(3,:) = [-0.5774 -0.5774 1 1];
    gauss_points2(4,:) = [-0.5774 0.5574 1 1];
    gauss_points2(5,:) = [0.5774 -0.5774 1 1];
    gauss_points2(6,:) = [0.5774 0.5774 1  1];
    gauss_points2(7,:) = [0.5774 -0.5774 1 1];
    gauss_points2(8,:) = [0.5774 0.5774 1  1];
end



new_phis = zeros(8,4); 
phicount = 1;
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:8
                new_phis(phicount,l) = (1+(-1)^(i)*gauss_points2(l,1))/2*...
                (1+(-1)^(k)*gauss_points2(l,2))/2*(1+(-1)^(j)*...
                gauss_points2(l,3))/2;
            end
            phicount = phicount + 1;
        end
    end
end




new_dphidxi = zeros(8,4);
new_dphidzeta = zeros(8,4);
new_dphideta = zeros(8,4);
xicount = 1;
for i =1:2
    for j=1:2
        for k=1:2
            for l=1:8
                new_dphidxi(xicount,l) = (-1)^(i)*0.5*(1+(-1)^(k)...
                *gauss_points2(l,2))*0.5*(1+(-1)^(j)*gauss_points2(l,3))*0.5; 
                
            end
            xicount = xicount + 1;
        end
    end
end
% Make a table of dphi/deta

etacount = 1;
for i=1:2
    for j=1:2
        for k = 1:2
            for l=1:8
                new_dphideta(etacount,l) = 0.5*(-1)^(k)*...
                    (1+(-1)^(i)*gauss_points2(l,1))*0.5*(1+(-1)^(j)*...
                    gauss_points2(l,3))*0.5;
                
            end
            etacount= etacount+1;
        end
    end
end
% Create a table of dphi/dzeta
zetacount = 1;
for i=1:2
    for j=1:2
        for k=1:2
            for l =1:8
                new_dphidzeta(zetacount, l) = 0.5*(-1)^(j)*(1+(-1)^(i)*...
                    gauss_points2(l,1))*0.5*(1+(-1)^(k)*gauss_points2(l,2))*...
                    0.5;
            end
            zetacount = zetacount + 1;
        end
    end
end

% Get the array of global node numbers per element (elx, ely, elz)
NOP_array = NOP(nny, nnz, nnx);
% Implement a Neumann boundary by adding the boundary integral to the
% b vector (passed in to this function as bvalues

% This applies the Neumann boundary condition to the entire boundary. USing
% the way the convection-diffusion code is structured now, we just want to
% do that on the particular element we have now
element = 1;
for elx=1:nnx-1
    for ely=1:nny-1
        for elz=1:nnz-1
            % get node coords for this elemnet
            nodes = NOP_array(element,:);
            % Find out if this element is on the boundary of interest
            if  onBoundary(direction, boundary_coordinates, n_coords(nodes,:));
                % We also have to extract the appropriate gp - not just
                % 1,2,3,4
                gps = [1 2 3 4 5 6 7 8];
                boundaries = boundary_node_matrix(nodes,:);
                gps = gps(boundaries(:,boundary));
                % Now, extract the nodes that are on the boundary
                % First get the boundary matrix
                
                real_nodes = nodes(boundaries(:,boundary));
                real_phis = new_phis(boundaries(:,boundary),:);
              for gp=1:4 % only 4 gps for a linear, quad, 2D element
                
                
                
                gp_weight = gauss_points2(gp,4); % column for is weight of gp
                % The above line must be changed if we switch to non
                % bilinear basis functions
                phi = real_phis(:,gps(gp)); 
                
                % Determine whether we need dphidxi, dphideta, dphidzeta
                if boundary ==1 || boundary ==2
                    dphidxis = new_dphidxi(boundaries(:,boundary),gps(gp)); % get every value of dphidxi at this gp
                 % This is a column vector
                    
                
                    dphisdzetas = new_dphidzeta(boundaries(:,boundary),gps(gp));
                    
                    % Need the Jacobian of transformation to evaluate integral
                    % in the computational space
                    dxdxi= sum(n_coords(real_nodes,1).*dphidxis);
                    dxdzeta = sum(n_coords(real_nodes,1).*dphisdzetas);
                    dzdxi = sum(n_coords(real_nodes,3).*dphidxis);
                    dzdzeta = sum(n_coords(real_nodes,3).*dphisdzetas);
                    Jacobian = [dxdxi dxdzeta; dzdxi dzdzeta];
                elseif boundary ==3 || boundary == 4
                     dphisdzetas = new_dphidzeta(boundaries(:,boundary),gps(gp));
                     dphisdetas = new_dphideta(boundaries(:,boundary),gps(gp));
                     dydeta = sum(n_coords(real_nodes,2).*dphisdetas);
                     dydzeta = sum(n_coords(real_nodes,2).*dphisdzetas);
                     dzdeta = sum(n_coords(real_nodes,3).*dphisdetas);
                     dzdzeta = sum(n_coords(real_nodes,3).*dphisdzetas);
                     Jacobian = [dydeta dydzeta; dzdeta dzdzeta];
                else
                     dphidxis = new_dphidxi(boundaries(:,boundary),gps(gp));
                     dphisdetas = new_dphideta(boundaries(:,boundary),gps(gp));
                     dxdxi= sum(n_coords(real_nodes,1).*dphidxis);
                     dxdeta = sum(n_coords(real_nodes,1).*dphisdetas);
                     dydxi = sum(n_coords(real_nodes,2).*dphidxis);
                     dydeta = sum(n_coords(real_nodes,2).*dphisdetas);
                     Jacobian = [dxdxi dxdeta; dydxi dydeta];
                end
                J = det(Jacobian);
                
                
                    % Get all the phis that we iterate over 
                        for node=1:4 % 4 is number of local nodes
                    % 
                            
                        
                       
                         bvalues(real_nodes(node)) = ...
                            bvalues(real_nodes(node))  +...
                            gp_weight*J*(value)*phi(node);
                            
                       
                        end
                        
                       
              end
            end
            element = element + 1;
        end
    end
end

b = bvalues;
end