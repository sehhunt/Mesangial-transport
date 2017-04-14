function [ plasma_flux ] = calculate_plasma_flux_dir( Pvector, ...
    km, mu, bounds, dir, meshx, meshy, meshz  ) %#codegen
%calculate_plasma_flux_dir Returns the plasma flux from any coordinate
%aligned surface, calculated at cell centers
%   Pvector is the vector of pressure at each node, km is the mesangial 
%   matrix permeability and mu is the fluid viscosity. bounds is the value
%   of x,y, or z at the boundary we want. dir is a vector with 1 (high   
%   boundary) or -1 in the position of an (x,y,z) vector that we want.
%   Positive flux means moving into the capillary of interest, negative
%   flux means moving out of the capillary of interest (loss of fluid).

% Find out which boundary we are checking

center_point = [dir 1]; % We are caluclating the flux at center of the
% surface 
    


num_phis = 8;
%num_gp = size(center_point,1);

dphisdxis = get_dphidxi(center_point,num_phis);
dphisdetas = get_dphideta(center_point, num_phis);
dphisdzetas = get_dphidzeta(center_point, num_phis);
node_coord_matrix = get_node_coords(meshx, meshy, meshz);

NYNodes = size(meshy,1);
NXNodes = size(meshx,2);
NZNodes = size(meshz,3);

dx = meshx(1,2,1)-meshx(1,1,1);
dy = meshy(2,1,1)-meshy(1,1,1);
dz = meshz(1,1,2)-meshz(1,1,1);

NOP_array = NOP(NYNodes, NZNodes, NXNodes);
element = 1;

if isequal(abs(dir(3)),1)
    plasma_flux = zeros(NYNodes-1,1);
else 
    plasma_flux = zeros(NZNodes-1,1);
end

for elx = 1:NXNodes-1
    for ely = 1:NYNodes-1
       for elz = 1:NZNodes-1
          el_nodes = NOP_array(element,:);
 
          if isequal(dir(1),-1)
              if onBoundary([true false false],bounds,node_coord_matrix(el_nodes,:))
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
                  area_conv = [dydeta dydzeta; dzdeta dzdzeta];
                  J = det(area_conv);
                  dphi = zeros(num_phis,3);
                  for i =1:num_phis
                      dphi(i,:) = conv_matrix\[dphisdxis(i); dphisdetas(i);...
                          dphisdzetas(i)];
                  end
                  
                  dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                                    
                  plasma_flux(elz) = plasma_flux(elz) + ...
                      km/mu*dpdx*abs(J)*4; %deta = 2 dzeta = 2
              end
          elseif isequal(dir(1),1)
              if onBoundary([true false false],bounds,node_coord_matrix(el_nodes,:))
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
                  
                  area_conv = [dydeta dydzeta; dzdeta dzdzeta];
                  J = det(area_conv);
                  
                  dphi = zeros(num_phis,3);
                  for i =1:num_phis
                      dphi(i,:) = conv_matrix\[dphisdxis(i); dphisdetas(i);...
                          dphisdzetas(i)];
                  end
                  
                  dpdx = sum(Pvector(el_nodes).*dphi(:,1));
                                    
                  plasma_flux(elz) = plasma_flux(elz) - ...
                      km/mu*dpdx*abs(J)*4;
              end
              
          elseif isequal(dir(2),-1)
              if onBoundary([false true false],bounds,node_coord_matrix(el_nodes,:))
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
                  
                  dpdy = sum(Pvector(el_nodes).*dphi(:,2));
                                    
                  plasma_flux(elz) = plasma_flux(elz) + ...
                      km/mu*dpdy*dx*dz;
              end
          elseif isequal(dir(2),1)
              if onBoundary([false true false],bounds,node_coord_matrix(el_nodes,:))
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
                  
                  dpdy = sum(Pvector(el_nodes).*dphi(:,2));
                                    
                  plasma_flux(elz) = plasma_flux(elz) - ...
                      km/mu*dpdy*dx*dz;
              end
          elseif isequal(dir(3),-1)
              if onBoundary([false false true],bounds,node_coord_matrix(el_nodes,:))
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
                  
                  dpdz = sum(Pvector(el_nodes).*dphi(:,3));
                                    
                  plasma_flux(ely) = plasma_flux(ely) + ...
                      km/mu*dpdz*dx*dy;
              end
          elseif isequal(dir(3),1)
              if onBoundary([false false true],bounds,node_coord_matrix(el_nodes,:))
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
                  
                  dpdz = sum(Pvector(el_nodes).*dphi(:,3));
                                    
                  plasma_flux(ely) = plasma_flux(ely) - ...
                      km/mu*dpdz*dx*dy;
              end
          end
          
          element = element + 1;
          
       end
    end
end

if isequal(abs(dir(1)),1)
    plasma_flux = plasma_flux./((NYNodes-1)*dy*dz);
elseif isequal(abs(dir(2)),1)
    plasma_flux = plasma_flux./(dz*dx*(NXNodes-1));
else
    plasma_flux = plasma_flux./(dx*dy*(NXNodes-1));
end

end

