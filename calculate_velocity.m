function [ velocity, node_coords ] = calculate_velocity( Pvector,x,y,z,...
    ybounds, km,kbm,mu )
% This function calculates the velocity based on the pressure field 
% (Pvector), the x. y, and z coordinates (formatted as outputs of 
% meshgrid), and the permeabilities of the mesangial matrix (km) and the
% basement membrane (kbm). Output is velocity, a matrix with column 1 the
% x component, column 2 the y component, and column 3 the z component of
% the velociy at node (row number). node_coords contains the x,y,z
% coordinates of all the gauss points, the points at which we find the
% velocities

% Instead of gauss points, use coordinates of element nodes in
% computational space
coord = 0.5774;
gauss_points = zeros(8,4);
gauss_points(1,:) = [-coord -coord, -coord  1];
gauss_points(2,:) = [-coord coord -coord 1];
gauss_points(3,:) = [-coord -coord coord 1];
gauss_points(4,:) = [-coord coord  coord 1];
gauss_points(5,:) = [coord -coord, -coord  1];
gauss_points(6,:) = [coord coord, -coord  1];
gauss_points(7,:) = [coord -coord, coord  1];
gauss_points(8,:) = [coord coord, coord  1];

dphisdxis = get_dphidxi(gauss_points, 8);
dphidzetas = get_dphidzeta(gauss_points, 8);
dphidetas = get_dphideta(gauss_points, 8);
phis = get_phis(gauss_points,8);
node_coord_matrix = get_node_coords(x,y,z);

num_ynodes = size(y,1);
num_xnodes = size(x,2);
num_znodes = size(z,3);
NOP_array = NOP(num_ynodes, num_znodes,num_xnodes);
element = 1;

velocity = zeros(length(Pvector),3);
node_coords = zeros(length(Pvector),3);
dphi = zeros(8,3);
for elx=1:num_xnodes-1
    for ely =1:num_ynodes-1
        for elz = 1:num_znodes-1
            
            el_nodes = NOP_array(element,:); % Returns the global node nums
            
            % Loop over each node, to calculate associated velocity
            for node=1:8
                               
                % Need to evaluate dxdxi, dydxi etc for this element
                dxdxi = sum(node_coord_matrix(el_nodes,1).*dphisdxis(:,node));
                dxdeta = sum(node_coord_matrix(el_nodes,1).*dphidetas(:,node));
                dxdzeta = sum(node_coord_matrix(el_nodes,1).*dphidzetas(:,node));
                
                dydxi = sum(node_coord_matrix(el_nodes,2).*dphisdxis(:,node));
                dydeta = sum(node_coord_matrix(el_nodes,2).*dphidetas(:,node));
                dydzeta = sum(node_coord_matrix(el_nodes,2).*dphidzetas(:,node));
                
                dzdxi = sum(node_coord_matrix(el_nodes,3).*dphisdxis(:,node));
                dzdeta = sum(node_coord_matrix(el_nodes,3).*dphidetas(:,node));
                dzdzeta = sum(node_coord_matrix(el_nodes,3).*dphidzetas(:,node));
                
                conversion_matrix = [dxdxi dydxi dzdxi; dxdeta dydeta dzdeta;...
                    dxdzeta dydzeta dzdzeta];
                
                for phi=1:8
                    dphi(phi,:) = conversion_matrix\[dphisdxis(phi,node); ...
                        dphidetas(phi,node); dphidzetas(phi,node)];
                end
                
                dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                dpdy = sum(Pvector(el_nodes).*dphi(:,2));
                dpdz = sum(Pvector(el_nodes).*dphi(:,3));
                
                % Now we need to check if the node is on the top y
                % boundary, b/c if it is, then we need to use kbm to
                % calculate the velocity, not km
                if ~onBoundary([false true false],ybounds(2),...
                        node_coord_matrix(el_nodes(node),:))
                    velocity(el_nodes(node),1) = ...
                        + -km/mu*dpdx;
                    velocity(el_nodes(node),2) = ...
                        -km/mu*dpdy;
                    velocity(el_nodes(node),3) = ...
                        - km/mu*dpdz;                    
                else
                    velocity(el_nodes(node),1) = ...
                        -kbm/mu*dpdx;
                    velocity(el_nodes(node),2) = ...
                        -kbm/mu*dpdy;
                    velocity(el_nodes(node),3) = ...
                        -kbm/mu*dpdz;
                    disp('true')
                end
                node_coords(el_nodes(node),1) = sum(...
                    node_coord_matrix(el_nodes,1).*phis(:,node));
                node_coords(el_nodes(node),2) = sum(...
                    node_coord_matrix(el_nodes,2).*phis(:,node));
                node_coords(el_nodes(node),3) = sum(...
                    node_coord_matrix(el_nodes,3).*phis(:,node));
            end
            element = element + 1;
        end
    end
end
end
