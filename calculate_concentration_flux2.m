function [ concentration_flux ] = calculate_concentration_flux2( Pvector, ...
    Cvector, km, mu, D, xbounds, islow, ...
    meshx, meshy, meshz )
%calculate_concentration_flux Calculates the total flux (diffusive +
%convective) out of the x boundary surface
%   Pvector is the vector of pressure at each node, Cvector is the vector
%   of concentration at each node. km is the mesangial matrix permeability
%   and mu is the fluid viscosity. xbounds is the value of x at the
%   boundary we want. islow is 1 (true) if we want the low x boundary, and
%   0 (false) if we want the high x boundary. direction is +1 if the flow
%   is in the positive z-direction. Direction is -1 if the flow is the
%   negative z-direction.
%   h: width of mesangial channel
%   Lz: length of capillary (z-direction)
%   D: diffusion coefficient for molecule of interest


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
phi = get_phis(center_point, num_phis);
node_coord_matrix = get_node_coords(meshx, meshy, meshz);

NYNodes = size(meshy,1);
NXNodes = size(meshx,2);
NZNodes = size(meshz,3);

NOP_array = NOP(NYNodes, NZNodes, NXNodes);
yspacing = meshy(2,1,1)-meshy(1,1,1);
dz = meshz(1,1,2)-meshz(1,1,1);
element = 1;

concentration_flux = zeros(NZNodes-1,1);

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
                  dcdx = sum(Cvector(el_nodes).*dphi(:,1));
                  center_c = sum(Cvector(el_nodes).*phi);
                  
                  concentration_flux(elz) = concentration_flux(elz) - ...
                      yspacing*dz*(D*dcdx + km/mu*dpdx*center_c);
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
                  dcdx = sum(Cvector(el_nodes).*dphi(:,1));
                  center_c = sum(Cvector(el_nodes).*phi);
                  
                  concentration_flux(elz) = concentration_flux(elz) +...
                      yspacing*dz*( ...
                      D*dcdx + km/mu*dpdx*center_c);
              end

          end
          
          element = element + 1;
          
       end
    end
end

concentration_flux = concentration_flux./(dz*yspacing*(NYNodes-1));
end

