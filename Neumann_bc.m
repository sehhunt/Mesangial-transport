% Apply a Neumann BC that depends on the dependent variable to boundary 2.
% A - stiffness matrix from Darcy problem (output of calculate_pressure) b:
% boundary condition vector, Cvector - stiffness matrix from diffusion
% problem (output of calculate_pressure), kbm - basement membrane matrix
% permeability, bm width - width of basement membrane, boundary_cond -
% pressure boundary condition on boundary 2 (y=H)

function [Amatrix, bvector] = Neumann_bc(A,b,Cvector,num_nodesx, ...
    num_nodesy, num_nodesz, mu, kbm,bm_width, boundary_cond, ybounds, ...
    xbounds, n_coords, varargin)
    % Construct boundary gauss points
    gauss_points2(1,:) = [-0.5774 1 -0.5774 1];
    gauss_points2(2,:) = [-0.5774 1 -0.5774 1];
    gauss_points2(3,:) = [-0.5774 1 0.5774 1];
    gauss_points2(4,:) = [-0.5774 1 0.5774 1];
    gauss_points2(5,:) = [0.5774 1 -0.5774 1];
    gauss_points2(6,:) = [0.5774 1 -0.5774 1];
    gauss_points2(7,:) = [0.5774 1 0.5774 1];
    gauss_points2(8,:) = [0.5774 1 0.5774 1];
    
    num_phis = 8;
    dphidxi = get_dphidxi(gauss_points2,num_phis);
    dphidzeta = get_dphidzeta(gauss_points2,num_phis);
    phis = get_phis(gauss_points2,num_phis);
    
    if length(varargin) == 1
        issinglecap = varargin{1};
    else
        issinglecap = true;
    end
  % This loop goes through all the elements, determines if they are on
  % boundary 2. If not, it does nothing. If they are on boundary 2 then it
  % applies a Neumann boundary condition to those nodes.
  for elx=1:num_nodesx-1
      for ely=1:num_nodesy-1
          for elz =1:num_nodesz-1
              element_nodes = test_element_nodes(elx, ely, elz, num_nodesy, num_nodesz);
              % Failing here onBoundary returns false
              if issinglecap
                  if (onBoundary([false true false],ybounds(2),n_coords(element_nodes,:))...
                          && ~onBoundary([true false false], xbounds(1), n_coords(element_nodes,:)))
                      
                      for gp=2:2:8 % only 4 gauss points on a boundary (gps 2,4,6,8 for this one)
                          phis_we_need = [2 4 6 8];
                          dphidxis = dphidxi(phis_we_need,gp); % get every value of dphidxi at this gp
                          % This is a column vector
                          dphisdzetas = dphidzeta(phis_we_need,gp);
                          real_phis_at_gp = phis(phis_we_need,gp);
                          % Calculate components of Jacobian
                          real_nodes = element_nodes(phis_we_need);
                          dxdxi= sum(n_coords(real_nodes,1).*dphidxis);
                          dxdzeta = sum(n_coords(real_nodes,1).*dphisdzetas);
                          dzdxi = sum(n_coords(real_nodes,3).*dphidxis);
                          dzdzeta = sum(n_coords(real_nodes,3).*dphisdzetas);
                          Jacobian = [dxdxi dxdzeta; dzdxi dzdzeta];
                          J = det(Jacobian);
                          
                          for i=1:4
                              for j=1:4
                                  A(real_nodes(i), real_nodes(j)) = ...
                                      A(real_nodes(i),real_nodes(j)) ...
                                      - (real_phis_at_gp(i)*real_phis_at_gp(j))*abs(J)...
                                      *kbm/(mu*bm_width);
                                  % CHANGE that needs to be bm_width no mes_width
                              end
                              bc_function = boundary_cond{2};
                              
                              b(real_nodes(i)) = b(real_nodes(i)) - abs(J)*kbm...
                                  /(mu*bm_width)*real_phis_at_gp(i)*...
                                  (bc_function(Cvector(real_nodes(i))));
                              
                              %end
                          end
                          
                      end
                

                  end
              else
                  if onBoundary([false true false],ybounds(2),n_coords(element_nodes,:))...
                          && ~onBoundary([true false false], xbounds(1), n_coords(element_nodes,:))...
                          ...
                          && ~onBoundary([true false false], xbounds(2), ...
                          n_coords(element_nodes,:))
                      
                      for gp=2:2:8 % only 4 gauss points on a boundary (gps 2,4,6,8 for this one)
                          phis_we_need = [2 4 6 8];
                          dphidxis = dphidxi(phis_we_need,gp); % get every value of dphidxi at this gp
                          % This is a column vector
                          dphisdzetas = dphidzeta(phis_we_need,gp);
                          real_phis_at_gp = phis(phis_we_need,gp);
                          % Calculate components of Jacobian
                          real_nodes = element_nodes(phis_we_need);
                          dxdxi= sum(n_coords(real_nodes,1).*dphidxis);
                          dxdzeta = sum(n_coords(real_nodes,1).*dphisdzetas);
                          dzdxi = sum(n_coords(real_nodes,3).*dphidxis);
                          dzdzeta = sum(n_coords(real_nodes,3).*dphisdzetas);
                          Jacobian = [dxdxi dxdzeta; dzdxi dzdzeta];
                          J = det(Jacobian);
                          
                          for i=1:4
                              for j=1:4
                                  A(real_nodes(i), real_nodes(j)) = ...
                                      A(real_nodes(i),real_nodes(j)) ...
                                      - (real_phis_at_gp(i)*real_phis_at_gp(j))*abs(J)...
                                      *kbm/(mu*bm_width);
                                  % CHANGE that needs to be bm_width no mes_width
                              end
                              bc_function = boundary_cond{2};
                              
                              b(real_nodes(i)) = b(real_nodes(i)) - abs(J)*kbm...
                                  /(mu*bm_width)*real_phis_at_gp(i)*...
                                  (bc_function(Cvector(real_nodes(i))));
                              
                              %end
                          end
                          
                      end
                

                  end
              end
          end
      end
      
  end
Amatrix = A;
bvector = b;
end