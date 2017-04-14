function [ plasma_flux ] = calculate_plasma_flux( Pvector, ...
    km, mu, xbounds, islow, meshx, meshy, meshz  ) %#codegen
%calculate_plasma_flux Returns the plasma flux on an x-surface, calculated
%at cell centers
%   Pvector is the vector of pressure at each node, km is the mesangial 
%   matrix permeability and mu is the fluid viscosity. xbounds is the value
%   of x at the boundary we want. islow is 1 (true) if we want the low x 
%   boundary, and 0 (false) if we want the high x boundary. 
%   Positive flux means moving into the capillary of interest, negative
%   flux means moving out of the capillary of interest (loss of fluid).

% Find out which boundary we are checking
if islow
    center_point = [-1 0 0 1]; % We are caluclating the flux at center of the
  % surface (x boundary, low or high)
else
    center_point = [1 0 0 1];
end

num_phis = 8;
%num_gp = size(center_point,1);

dphisdxis = get_dphidxi(center_point,num_phis);
dphisdetas = get_dphideta(center_point, num_phis);
dphisdzetas = get_dphidzeta(center_point, num_phis);
node_coord_matrix = get_node_coords(meshx, meshy, meshz);

NYNodes = size(meshy,1);
NXNodes = size(meshx,2);
NZNodes = size(meshz,3);

NOP_array = NOP(NYNodes, NZNodes, NXNodes);
element = 1;

plasma_flux = zeros(NZNodes-1,1);

for elx = 1:NXNodes-1
    for ely = 1:NYNodes-1
       for elz = 1:NZNodes-1
          el_nodes = NOP_array(element,:);
 
          if islow
              if onBoundary([true false false],xbounds,node_coord_matrix(el_nodes,:))
                  dxdxi = sum(node_coord_matrix(el_nodes,1).*dphisdxis);
                  dxdeta = sum(node_coord_matrix(el_nodes,1).*dphisdetas);
                  dxdzeta = sum(node_coord_matrix(el_nodes,1).*dphisdzetas);
                  dydxi = sum(node_coord_matrix(el_nodes,2).*dphisdxis);
                  dydeta = sum(node_coord_matrix(el_nodes,2).*dphisdetas);
                  dydzeta = sum(node_coord_matrix(el_nodes,2).*dphisdzetas);
                  dzdxi = sum(node_coord_matrix(el_nodes,3).*dphisdxis);
                  dzdeta = sum(node_coord_matrix(el_nodes,3).*dphisdetas);
                  dzdzeta = sum(node_coord_matrix(el_nodes,3).*dphisdzetas);
                  
                  conv_matrix = [dxdxi dydxi dzdxi; dxdeta dydeta dzdeta;...
                      dxdzeta dydzeta dzdzeta];
                  
                  dphi = zeros(num_phis,3);
                  for i =1:num_phis
                      dphi(i,:) = conv_matrix\[dphisdxis(i); dphisdetas(i);...
                          dphisdzetas(i)];
                  end
                  
                  dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                                    
                  plasma_flux(elz) = plasma_flux(elz) - ...
                      km/mu*dpdx;
              end
          else
              if onBoundary([true false false],xbounds,node_coord_matrix(el_nodes,:))
                  dxdxi = sum(node_coord_matrix(el_nodes,1).*dphisdxis);
                  dxdeta = sum(node_coord_matrix(el_nodes,1).*dphisdetas);
                  dxdzeta = sum(node_coord_matrix(el_nodes,1).*dphisdzetas);
                  dydxi = sum(node_coord_matrix(el_nodes,2).*dphisdxis);
                  dydeta = sum(node_coord_matrix(el_nodes,2).*dphisdetas);
                  dydzeta = sum(node_coord_matrix(el_nodes,2).*dphisdzetas);
                  dzdxi = sum(node_coord_matrix(el_nodes,3).*dphisdxis);
                  dzdeta = sum(node_coord_matrix(el_nodes,3).*dphisdetas);
                  dzdzeta = sum(node_coord_matrix(el_nodes,3).*dphisdzetas);
                  
                  conv_matrix = [dxdxi dydxi dzdxi; dxdeta dydeta dzdeta;...
                      dxdzeta dydzeta dzdzeta];
                  
                  dphi = zeros(num_phis,3);
                  for i =1:num_phis
                      dphi(i,:) = conv_matrix\[dphisdxis(i); dphisdetas(i);...
                          dphisdzetas(i)];
                  end
                  
                  dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                                    
                  plasma_flux(elz) = plasma_flux(elz) + ...
                      km/mu*dpdx;
              end

          end
          
          element = element + 1;
          
       end
    end
end

plasma_flux = plasma_flux./(NYNodes-1);
end

