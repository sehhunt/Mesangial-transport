% This function takes in a matrix (which is the output of a finite element
% code), and the mesh of the domain for that matrix, and calculates flux
% over the various surfaces of that domain

function [Integral,bxlow,bxhigh,bylow,byhigh,bzlow,bzhigh,non_zero_lowy] = ...
    calculate_flux(Pvector, k,mu, x,y,z, ...
    xbounds,ybounds,zbounds)
% Create gauss points
gauss_points = zeros(8,4);
gauss_points(1,:) = [-0.5774 -0.5774, -0.5774  1];
gauss_points(2,:) = [-0.5774 0.5774 -0.5774 1];
gauss_points(3,:) = [-0.5774 -0.5774 0.5774 1];
gauss_points(4,:) = [-0.5774 0.5774  0.5774 1];
gauss_points(5,:) = [0.5774 -0.5774, -0.5774  1];
gauss_points(6,:) = [0.5774 0.5774, -0.5774  1];
gauss_points(7,:) = [0.5774 -0.5774, 0.5774  1];
gauss_points(8,:) = [0.5774 0.5774, 0.5774  1];

%phis = get_phis(gauss_points, 8);
dphisdxis = get_dphidxi(gauss_points, 8);
dphidzetas = get_dphidzeta(gauss_points, 8);
dphidetas = get_dphideta(gauss_points, 8);
node_coord_matrix = get_node_coords(x,y,z);

lowxsurface_gauss_points = zeros(4,4);
lowxsurface_gauss_points(1,:) = [-1 -0.5774 -0.5774 1];
lowxsurface_gauss_points(2,:) = [-1 0.5774 -0.5774 1];
lowxsurface_gauss_points(3,:) = [-1 -0.5774 0.5774 1];
lowxsurface_gauss_points(4,:) = [-1 0.5774 0.5774 1];

% Get the derivatives on the low x surface
lowx_dphideta = get_dphideta(lowxsurface_gauss_points,8);
lowx_dphidzeta = get_dphidzeta(lowxsurface_gauss_points,8);
lowx_dphidxi = get_dphidxi(lowxsurface_gauss_points,8);

highxsurface_gauss_points = zeros(4,4);
highxsurface_gauss_points(1,:) = [1 -0.5774 -0.5774 1];
highxsurface_gauss_points(2,:) = [1 0.5774 -0.5774 1];
highxsurface_gauss_points(3,:) = [1 -0.5774 0.5774 1];
highxsurface_gauss_points(4,:) = [1 0.5774 0.5774 1];

highx_dphideta = get_dphideta(highxsurface_gauss_points,8);
highx_dphidzeta = get_dphidzeta(highxsurface_gauss_points,8);
highx_dphidxi = get_dphidxi(highxsurface_gauss_points,8);

lowysurface_gauss_points = zeros(4,4);

lowysurface_gauss_points(1,:) = [-0.5774 -1 -0.5774 1];
lowysurface_gauss_points(2,:) = [-0.5774 -1 0.5774 1];
lowysurface_gauss_points(3,:) = [0.5774 -1 -0.5774 1];
lowysurface_gauss_points(4,:) = [0.5774 -1 0.5774 1];

lowy_dphideta = get_dphideta(lowysurface_gauss_points,8);
lowy_dphidzeta = get_dphidzeta(lowysurface_gauss_points,8);
lowy_dphidxi = get_dphidxi(lowysurface_gauss_points,8);

highysurface_gauss_points(1,:) = [-0.5774 1 -0.5774 1];
highysurface_gauss_points(2,:) = [-0.5774 1 0.5774 1];
highysurface_gauss_points(3,:) = [0.5774 1 -0.5774 1];
highysurface_gauss_points(4,:) = [0.5774 1 0.5774 1];

highy_dphideta = get_dphideta(highysurface_gauss_points,8);
highy_dphidzeta = get_dphidzeta(highysurface_gauss_points,8);
highy_dphidxi = get_dphidxi(highysurface_gauss_points,8);

lowz_gauss_points(1,:) = [-0.5774 -0.5774 -1 1];
lowz_gauss_points(2,:) = [-0.5774 0.5774 -1 1];
lowz_gauss_points(3,:) = [0.5774 -0.5774 -1 1];
lowz_gauss_points(4,:) = [0.5774 0.5774 -1 1];

lowz_dphideta = get_dphideta(lowz_gauss_points,8);
lowz_dphidzeta = get_dphidzeta(lowz_gauss_points,8);
lowz_dphidxi = get_dphidxi(lowz_gauss_points,8);

highz_gauss_points(1,:) = [-0.5774 -0.5774 1 1];
highz_gauss_points(2,:) = [-0.5774 0.5774 1 1];
highz_gauss_points(3,:) = [0.5774 -0.5774 1 1];
highz_gauss_points(4,:) = [0.5774 0.5774 1 1];

highz_dphideta = get_dphideta(highz_gauss_points,8);
highz_dphidzeta = get_dphidzeta(highz_gauss_points,8);
highz_dphidxi = get_dphidxi(highz_gauss_points,8);

NNODES = length(node_coord_matrix);
boundaries = false (NNODES, 6);
for i=1:NNODES
    boundaries(i,:) = [onBoundary([false true false], ybounds(1), ...
        node_coord_matrix(i,:)) onBoundary([false true false], ybounds(2), ...
        node_coord_matrix(i,:)) onBoundary([true false false], xbounds(1), ...
        node_coord_matrix(i,:)) onBoundary([true false false], xbounds(2), ...
        node_coord_matrix(i,:)) onBoundary([false false true], zbounds(1), ...
        node_coord_matrix(i,:)) onBoundary([false false true], zbounds(2), ...
        node_coord_matrix(i,:)) ];
end
Integral = 0;
num_ynodes = size(y,1);
num_xnodes = size(x,2);
num_znodes = size(z,3);
NOP_array = NOP(num_ynodes, num_znodes,num_xnodes);
element = 1;
bxlow = 0;
bxhigh = 0;
bylow = 0;
byhigh = 0;
bzlow = 0;
bzhigh = 0;
nz_lowy_count = 1;
for elx=1:num_xnodes-1
    for ely=1:num_ynodes-1
        for elz=1:num_znodes-1
            % Need to do separate calculations for each boundary
            el_nodes = NOP_array(element,:); % Returns global node numbers 
                % of the nodes in the current element
                % Test to see if on low_x boundary
                
                % Per element, there should be one value for del P. Use
                % this on all appropriate boundaries
                if onBoundary([true false false],xbounds(1),...
                        node_coord_matrix(el_nodes,:))
                    lowx_phis = [1 2 3 4]; % Only these nodes are on this boundary
                    % I need to evaluate both dp/dx and dydz on this
                    % boundary. In order to get dp/dx I need to use all of
                    % the dphis - but where are these dphis evaluated?
                    for gp = 1:4
                        
                        real_nodes = el_nodes(lowx_phis);
%                         real_dphisdetas = dphidetas(lowx_phis,gp);
%                         real_dphisdzetas = dphidzetas(lowx_phis,gp);
%                         real_dphisdxis = dphisdxis(lowx_phis,gp);
                        dydeta = sum(node_coord_matrix(el_nodes,2).*...
                            lowx_dphideta(:,gp));
                        dydzeta = sum(node_coord_matrix(el_nodes,2).*...
                            lowx_dphidzeta(:,gp));
                        dzdeta = sum(node_coord_matrix(el_nodes,3).*...
                            lowx_dphideta(:,gp));
                        dzdzeta = sum(node_coord_matrix(el_nodes,3).*...
                            lowx_dphidzeta(:,gp));
                        dydxi = sum(node_coord_matrix(el_nodes,2).*...
                            lowx_dphidxi(:,gp));
                        dzdxi = sum(node_coord_matrix(el_nodes,3).*...
                            lowx_dphidxi(:,gp));
                        dxdxi = sum(node_coord_matrix(el_nodes,1).*...
                            lowx_dphidxi(:,gp));
                        dxdeta = sum(node_coord_matrix(el_nodes,1).*...
                            lowx_dphideta(:,gp));
                        dxdzeta = sum(node_coord_matrix(el_nodes,1).*...
                            dphidzetas(:,gp));
                        
                        Jacobian2 = [dxdxi dydxi dzdxi; dxdeta dydeta dzdeta;...
                            dxdzeta dydzeta dzdzeta];
                        
                        dphi = zeros(8,3);
                        for phi = 1:8
                            dphi(phi,:) = Jacobian2\[lowx_dphidxi(phi,gp);...
                                lowx_dphideta(phi,gp); lowx_dphidzeta(phi,gp)];
                        end
                        dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                        dydeta_surf = sum(node_coord_matrix(real_nodes,2).*...
                            lowx_dphideta(lowx_phis,gp));
                        dydzeta_surf = sum(node_coord_matrix(real_nodes,2).*...
                            lowx_dphidzeta(lowx_phis,gp));
                        dzdeta_surf = sum(node_coord_matrix(real_nodes,3).*...
                            lowx_dphideta(lowx_phis,gp));
                        dzdzeta_surf = sum(node_coord_matrix(real_nodes,3).*...
                            lowx_dphidzeta(lowx_phis,gp));
                        Jacobian = [dydeta_surf dydzeta_surf; ...
                            dzdeta_surf dzdzeta_surf];
                        J = det(Jacobian);
                        
                            bxlow = bxlow  ...
                                + J*k/mu*dpdx;
                       
                    end
                 
                    % Otherwise, check if it's on the xhigh boundary
                elseif onBoundary([true false false],xbounds(2),...
                        node_coord_matrix(el_nodes,:))
                    highx_phis = [5 6 7 8];
                    
                    for gp = 1:4
                         real_nodes = el_nodes(highx_phis);
%                         real_dphisdetas = dphidetas(highx_phis,gp);
%                         real_dphisdzetas = dphidzetas(highx_phis,gp);
%                         real_dphisdxis = dphisdxis(highx_phis,gp);
                        dydeta = sum(node_coord_matrix(el_nodes,2).*...
                            highx_dphideta(:,gp));
                        dydzeta = sum(node_coord_matrix(el_nodes,2).*...
                            highx_dphidzeta(:,gp));
                        dzdeta = sum(node_coord_matrix(el_nodes,3).*...
                            highx_dphideta(:,gp));
                        dzdzeta = sum(node_coord_matrix(el_nodes,3).*...
                            highx_dphidzeta(:,gp));
                        dydxi = sum(node_coord_matrix(el_nodes,2).*...
                            highx_dphidxi(:,gp));
                        dzdxi = sum(node_coord_matrix(el_nodes,3).*...
                            highx_dphidxi(:,gp));
                        dxdxi = sum(node_coord_matrix(el_nodes,1).*...
                            highx_dphidxi(:,gp));
                        dxdeta = sum(node_coord_matrix(el_nodes,1).*...
                            highx_dphideta(:,gp));
                        dxdzeta = sum(node_coord_matrix(el_nodes,1).*...
                            highx_dphidzeta(:,gp));
                        
                        Jacobian2 = [dxdxi dxdeta dxdzeta; dydxi dydeta dydzeta;...
                            dzdxi dzdeta dzdzeta];
                        
                        dphi = zeros(8,3);
                        for phi = 1:8
                            dphi(phi,:) = Jacobian2\[highx_dphidxi(phi,gp); ...
                                highx_dphideta(phi,gp); ...
                                highx_dphidzeta(phi,gp)];
                        end
                        dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                        dydeta_surf = sum(node_coord_matrix(real_nodes,2).*...
                            highx_dphideta(highx_phis,gp));
                        dydzeta_surf = sum(node_coord_matrix(real_nodes,2).*...
                            highx_dphidzeta(highx_phis,gp));
                        dzdeta_surf = sum(node_coord_matrix(real_nodes,3).*...
                            highx_dphideta(highx_phis,gp));
                        dzdzeta_surf = sum(node_coord_matrix(real_nodes,3).*...
                            highx_dphidzeta(highx_phis,gp));
                        Jacobian = [dydeta_surf dydzeta_surf; ...
                            dzdeta_surf dzdzeta_surf];
                        J = det(Jacobian);
                        
                            bxhigh = bxhigh + ...
                                -J*k/mu*dpdx;
                        
                    end
                end
                % Check to see if element is on the ylow boundary
                if onBoundary([false true false],ybounds(1),...
                        node_coord_matrix(el_nodes,:))
                    lowy_phis = [1 3 5 7];
                    real_nodes = el_nodes(lowy_phis);
                    for gp = 1:4
%                         real_dphisdxis = dphisdxis(lowy_phis,gp);
%                         real_dphidzetas = dphidzetas(lowy_phis,gp);
%                         real_dphidetas = dphidetas(lowy_phis,gp);
                        dxdxi = sum(node_coord_matrix(el_nodes,1).*...
                            lowy_dphidxi(:,gp));
                        dxdzeta = sum(node_coord_matrix(el_nodes,1).*...
                            lowy_dphidzeta(:,gp));
                        dzdxi = sum(node_coord_matrix(el_nodes,3).*...
                            lowy_dphidxi(:,gp));
                        dzdzeta = sum(node_coord_matrix(el_nodes,3).*...
                            lowy_dphidzeta(:,gp));
                        dxdeta = sum(node_coord_matrix(el_nodes,1).*...
                            lowy_dphideta(:,gp));
                        dzdeta = sum(node_coord_matrix(el_nodes,3).*...
                            lowy_dphideta(:,gp));
                        dydxi = sum(node_coord_matrix(el_nodes,2).*...
                            lowy_dphidxi(:,gp));
                        dydeta = sum(node_coord_matrix(el_nodes,2).*...
                            lowy_dphideta(:,gp));
                        dydzeta = sum(node_coord_matrix(el_nodes,2).*...
                            lowy_dphidzeta(:,gp));
                       
                        Jacobian2 = [dxdxi dxdeta dxdzeta; dydxi dydeta dydzeta;...
                            dzdxi dzdeta dzdzeta];
                       
                        dphi = zeros(8,3);
                        for phi = 1:8
                            dphi(phi,:) = Jacobian2\[lowy_dphidxi(phi,gp); ...
                                lowy_dphideta(phi, gp);...
                                lowy_dphidzeta(phi, gp)];
                        end
                        dpdy = sum(Pvector(el_nodes).*dphi(:,2));
                        if dpdy ~= 0
                            non_zero_lowy(nz_lowy_count,1) = dpdy;
                            non_zero_lowy(nz_lowy_count,2) = element;
                            nz_lowy_count = nz_lowy_count+1;
                        end
                        dxdxi_surf = sum(node_coord_matrix(real_nodes,1).*...
                            lowy_dphidxi(lowy_phis,gp));
                        dxdzeta_surf = sum(node_coord_matrix(real_nodes,1).*...
                            lowy_dphidzeta(lowy_phis,gp));
                        dzdxi_surf= sum(node_coord_matrix(real_nodes,3).*...
                            lowy_dphidxi(lowy_phis,gp));
                        dzdzeta_surf = sum(node_coord_matrix(real_nodes,3).*...
                            lowy_dphidzeta(lowy_phis,gp));
                        Jacobian = [dxdxi_surf dxdzeta_surf; ...
                            dzdxi_surf dzdzeta_surf];
                         J = det(Jacobian);
                        
                            bylow = bylow + ...
                                k/mu*J*dpdy;
                       
                    end
                % Check to see if element is on the yhigh boundary
                elseif onBoundary([false true false], ybounds(2),...
                        node_coord_matrix(el_nodes,:))
                    highy_phis = [2 4 6 8];
                    real_nodes = el_nodes(highy_phis);
                    for gp = 1:4
%                         real_dphisdxis = dphisdxis(highy_phis,gp);
%                         real_dphidzetas = dphidzetas(highy_phis,gp);
%                         real_dphidetas = dphidetas(lowy_phis,gp);
                        dxdxi = sum(node_coord_matrix(el_nodes,1).*...
                            highy_dphidxi(:,gp));
                        dxdzeta = sum(node_coord_matrix(el_nodes,1).*...
                            highy_dphidzeta(:,gp));
                        dzdxi = sum(node_coord_matrix(el_nodes,3).*...
                            highy_dphidxi(:,gp));
                        dzdzeta = sum(node_coord_matrix(el_nodes,3).*...
                            highy_dphidzeta(:,gp));
                        dxdeta = sum(node_coord_matrix(el_nodes,1).*...
                            highy_dphideta(:,gp));
                        dzdeta = sum(node_coord_matrix(el_nodes,3).*...
                            highy_dphideta(:,gp));
                        dydxi = sum(node_coord_matrix(el_nodes,2).*...
                            highy_dphidxi(:,gp));
                        dydeta = sum(node_coord_matrix(el_nodes,2).*...
                            highy_dphideta(:,gp));
                        dydzeta = sum(node_coord_matrix(el_nodes,2).*...
                            highy_dphidzeta(:,gp));
                        
                        Jacobian2 = [dxdxi dxdeta dxdzeta; dydxi dydeta dydzeta;...
                            dzdxi dzdeta dzdzeta];
                       
                        dphi = zeros(8,3);
                        for phi = 1:4
                            dphi(phi,:) = Jacobian2\[highy_dphidxi(phi,gp); ...
                                highy_dphideta(phi,gp); highy_dphidzeta(phi,gp)];
                        end
                        dpdy = sum(Pvector(el_nodes).*dphi(:,2));
                        dxdxi_surf = sum(node_coord_matrix(real_nodes,1).*...
                            highy_dphidxi(highy_phis,gp));
                        dxdzeta_surf = sum(node_coord_matrix(real_nodes,1).*...
                            highy_dphidzeta(highy_phis,gp));
                        dzdxi_surf = sum(node_coord_matrix(real_nodes,3).*...
                            highy_dphidxi(highy_phis,gp));
                        dzdzeta_surf = sum(node_coord_matrix(real_nodes,3).*...
                            highy_dphidzeta(highy_phis,gp));
                        Jacobian = [dxdxi_surf dxdzeta_surf; ...
                            dzdxi_surf dzdzeta_surf];
                         J = det(Jacobian);
                        
                            byhigh = byhigh - ...
                                k/mu*J*dpdy;
                        
                    end
                end
                % Check to see if the element is on the zlow boundary
                if onBoundary([false false true], zbounds(1),...
                    node_coord_matrix(el_nodes,:))
                    lowz_phis = [1 2 5 6];
                    real_nodes = el_nodes(lowz_phis);
                    for gp = 1:4
%                         real_dphisdxis = dphisdxis(lowz_phis,gp);
%                         real_dphidetas = dphidetas(lowz_phis,gp);
%                         real_dphidzetas = dphidzetas(lowz_phis,gp);
                        dydeta = sum(node_coord_matrix(el_nodes,2).*...
                            lowz_dphideta(:,gp));
                        dydxi = sum(node_coord_matrix(el_nodes,2).*...
                            lowz_dphidxi(:,gp));
                        dxdxi = sum(node_coord_matrix(el_nodes,1).*...
                            lowz_dphidxi(:,gp));
                        dxdeta = sum(node_coord_matrix(el_nodes,1).*...
                            lowz_dphideta(:,gp));
                        dydzeta = sum(node_coord_matrix(el_nodes,2).*...
                            lowz_dphidzeta(:,gp));
                        dxdzeta = sum(node_coord_matrix(el_nodes,1).*...
                            lowz_dphidzeta(:,gp));
                        dzdxi = sum(node_coord_matrix(el_nodes,3).*...
                            lowz_dphidxi(:,gp));
                        dzdeta = sum(node_coord_matrix(el_nodes,3).*...
                            lowz_dphideta(:,gp));
                        dzdzeta = sum(node_coord_matrix(el_nodes,3).*...
                            lowz_dphidzeta(:,gp));
                        
                        Jacobian2 = [dxdxi dxdeta dxdzeta;...
                            dydxi dydeta dydzeta; dzdxi dzdeta dzdzeta];
                       
                        dphi = zeros(8,3);
                        for phi = 1:8
                            dphi(phi,:) = Jacobian2\[lowz_dphidxi(phi,gp); ...
                                lowz_dphideta(phi,gp); lowz_dphidzeta(phi,gp)];
                        end
                        dpdz = sum(Pvector(el_nodes).*dphi(:,3));
                        dxdxi_surf = sum(node_coord_matrix(real_nodes,1).*...
                            lowz_dphidxi(lowz_phis,gp));
                        dxdeta_surf = sum(node_coord_matrix(real_nodes,1).*...
                            lowz_dphideta(lowz_phis,gp));
                        dydxi_surf = sum(node_coord_matrix(real_nodes,2).*...
                            lowz_dphidxi(lowz_phis,gp));
                        dydeta_surf = sum(node_coord_matrix(real_nodes,2).*...
                            lowz_dphidzeta(lowz_phis,gp));
                        Jacobian = [dxdxi_surf dxdeta_surf; ...
                            dydxi_surf dydeta_surf];
                        J = det(Jacobian);
                        
                            bzlow = bzlow + k/mu*J*dpdz;
                        
                    end
                    
                elseif onBoundary([false false true], zbounds(2),...
                        node_coord_matrix(el_nodes,:))
                    highz_phis = [3 4 7 8];
                    real_nodes = el_nodes(highz_phis);
                    for gp = 1:4
%                         real_dphisdxis = dphisdxis(highz_phis,gp);
%                         real_dphidetas = dphidetas(highz_phis,gp);
%                         real_dphidzetas = dphidzetas(highz_phis,gp);
                        dydeta = sum(node_coord_matrix(el_nodes,2).*...
                            highz_dphideta(:,gp));
                        dydxi = sum(node_coord_matrix(el_nodes,2).*...
                            highz_dphidxi(:,gp));
                        dxdxi = sum(node_coord_matrix(el_nodes,1).*...
                            highz_dphidxi(:,gp));
                        dxdeta = sum(node_coord_matrix(el_nodes,1).*...
                            highz_dphideta(:,gp));
                        dxdzeta = sum(node_coord_matrix(el_nodes,1).*...
                            highz_dphidzeta(:,gp));
                        dydzeta = sum(node_coord_matrix(el_nodes,2).*...
                            highz_dphidzeta(:,gp));
                        dzdxi = sum(node_coord_matrix(el_nodes,3).*...
                            highz_dphidxi(:,gp));
                        dzdeta = sum(node_coord_matrix(el_nodes,3).*...
                            highz_dphideta(:,gp));
                        dzdzeta = sum(node_coord_matrix(el_nodes,3).*...
                            highz_dphidzeta(:,gp));
                        
                        Jacobian2 = [dxdxi dxdeta dxdzeta;...
                            dydxi dydeta dydzeta; dzdxi dzdeta dzdzeta];
                        
                        dphi = zeros(8,3);
                        for phi = 1:8
                            dphi(phi,:) = Jacobian2\[highz_dphidxi(phi,gp); ...
                                highz_dphideta(phi,gp); highz_dphidzeta(phi,gp)];
                        end
                        dpdz = sum(Pvector(el_nodes).*dphi(:,3));
                        dxdxi_surf = sum(node_coord_matrix(real_nodes,1).*...
                            highz_dphidxi(highz_phis,gp));
                        dxdeta_surf = sum(node_coord_matrix(real_nodes,1).*...
                            highz_dphideta(highz_phis,gp));
                        dydxi_surf = sum(node_coord_matrix(real_nodes,2).*...
                            highz_dphidxi(highz_phis,gp));
                        dydeta_surf = sum(node_coord_matrix(real_nodes,2).*...
                            highz_dphideta(highz_phis,gp));
                        Jacobian = [dxdxi_surf dxdeta_surf; ...
                            dydxi_surf dydeta_surf];
                        J = det(Jacobian);
                        
                        
                            bzhigh = bzhigh - k/mu*J*dpdz;
                        
                    end
                end
                
                
                element = element + 1;
        end
    end
    
end

Integral = bxlow+bxhigh+bylow+byhigh+bzlow+bzhigh;
end