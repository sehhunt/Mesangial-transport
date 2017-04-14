% This code solves linked convection and diffusion problems
% It needs a mesh to solve over, and boundary conditions?
% xbounds, ybounds, and zbounds are vectors of length two giving the
% xvalues of the left, and right boundaries of the domain in xbounds, the
% bottom and top y vaules in ybounds (in that order), and the back and
% front z value (that order) in zbounds. vargargin takes the boundary
% condition functions, if they exist
function [P,C,CIg,loop_counter, varargout] = convection_diffusion(meshx, meshy, meshz, xbounds, ...
    ybounds, zbounds, P_bcs, C_bcs, Ig_bcs, bm_width, mu, kbm, km, Dalb,...
    DIgA, varargin) 
tic; % Time how long the code takes to run

% If meshx, meshy, and meshz were created using 3D meshgrid (which I am
% currently assuming, then they contain the x,y, and z coordinates of the
% discrete points of the grid. columns of meshy contain increasing
% coordinates, as rows are y value. rows of meshx contain increasing
% values, as coltumns are x coordinate. meshz increases in the third
% dimension (pages)
num_nodesx = size(meshx, 2);
num_nodesy = size(meshy, 1);
num_nodesz = size(meshz, 3);

new_meshx = meshx(:,:,1);
new_meshz = meshz(:,:,1);

for i =2:num_nodesz
    new_meshx = cat(1, new_meshx, meshx(:,:,i));
    new_meshz = cat(1, new_meshz, meshz(:,:,i));

end

meshx = new_meshx;
meshz = new_meshz;


% Table of x,y,z coordinates of each node, along with global number
% Layout: row number is global node number  x coord (col 1) y coord (col 2)
% z coord (col 3) 
node_coords = zeros(num_nodesx*num_nodesy*num_nodesz, 3);
NNODES = num_nodesx*num_nodesy*num_nodesz;
NELEMENTS = (num_nodesx-1)*(num_nodesy-1)*(num_nodesz-1);
for i=1:NNODES
    node_coords(i,1:3) = [meshx(i) meshy(i) meshz(i)];
end
% Matrix called boundaries with 6 columns that specify
% if the node is part of a boundary 1-6,
% boundary 1,2 are y=0,H, boundary 3,4 are x=0,L, and boundary 5,6 are
% z=0,W. If a node is on a boundary than the value of the matrix in that
% column is 1, otherwise it is zero

% This is the boundaries for the pressure term
boundaries = false(NNODES,6);

for i=1:NNODES
    boundaries(i,:) = [onBoundary([false true false], ybounds(1), ...
        node_coords(i,:)) onBoundary([false true false], ybounds(2), ...
        node_coords(i,:)) onBoundary([true false false], xbounds(1), ...
        node_coords(i,:)) onBoundary([true false false], xbounds(2), ...
        node_coords(i,:)) onBoundary([false false true], zbounds(1), ...
        node_coords(i,:)) onBoundary([false false true], zbounds(2), ...
        node_coords(i,:)) ];
end

if length(varargin) == 1
    issinglecap = varargin{1};
    
else
    issinglecap = true;
end
   
kappa = ones(NELEMENTS,1).*km;
                                                                                                                                                                                                                                                                                                                                                                                                                                                       
Amatrix = spalloc(NNODES, ...
    NNODES, 6*NNODES);

Cmatrix1 = spalloc(NNODES, NNODES, 6*NNODES);
Cmatrix2 = spalloc(NNODES, NNODES, 6*NNODES);
 Amatrix = calculate_pressure(Amatrix, num_nodesx,num_nodesy,num_nodesz,...
        node_coords,kappa,mu);
    Cvector = zeros(NNODES,1);
    
    % Construct the boundary condition vector
    b = zeros(NNODES, 1);
    NOP_array = NOP(num_nodesy,num_nodesz,num_nodesx);
    
    % Attempt to apply boundary conditions, take 2. This time we check for
    % the type of condition first - first applying neumann conditions, and
    % then applying Dirichlet conditions, which will override Neumann
    % conditions
        
    el = 1;
    for elx = 1:num_nodesx-1
        for ely = 1:num_nodesy-1
            for elz = 1:num_nodesz-1
                                
                % First check the low x boundary (usually x = 0)
                if onBoundary([true false false], xbounds(1),...
                        node_coords(NOP_array(el,:),:)) 
                    el_nodes = NOP_array(el,:);
                    % Is it a Drichlet condition? If so, then iterate by
                    % nodes. Otherwise, operate element-wise
                    bc = P_bcs{3};
                    
                    if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
%                     
                        b = one_el_neumann(3,xbounds,boundaries,node_coords,...
                                    bc{2},b,num_nodesy,num_nodesx,...
                                    num_nodesz,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        
                         b = one_el_neumann(3,xbounds,boundaries,node_coords,...
                                    bc{2},b,num_nodesy,num_nodesx,...
                                    num_nodesz,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                    
                elseif onBoundary([true false false], xbounds(2),...
                        node_coords(NOP_array(el,:),:)) 
                    % That just checked boundary 4
                    % Is there a Dirichlet condition on this boundary?
                    bc = P_bcs{4};
                    el_nodes = NOP_array(el,:);

                    if eq(bc{1},'N') && (~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)) 
                    
                        b = one_el_neumann(4,xbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(4,xbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                end
                
                % Now, check the y boundaries (1 & 2)
                 if onBoundary([false true false],ybounds(1),...
                         node_coords(NOP_array(el,:),:))
                     bc = P_bcs{1};
                        
                    if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
                        b = one_el_neumann(1,ybounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(1,ybounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                elseif onBoundary([false true false],ybounds(2),...
                        node_coords(NOP_array(el,:),:))
                     bc = P_bcs{2};
                     el_nodes = NOP_array(el,:);

                     if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
                          b = one_el_neumann(2,ybounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                       % Don't include this elseif, becasue the
                       % function-handle here is dealt with by Neumann_bc
%                      elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
%                          b = one_el_neumann(2,ybounds,boundaries,node_coords,...
%                             bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
%                             elx,ely,elz,Cvector,Amatrix...
%                                     ,kbm,bm_width,mu);
                     end
                end
                
                if onBoundary([false false true], zbounds(1),...
                        node_coords(NOP_array(el,:),:))
                    bc = P_bcs{5};
                    el_nodes = NOP_array(el,:);

                    if  eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0) 
                        b = one_el_neumann(5,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz...
                            ,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(5,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz...
                            ,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                elseif onBoundary([false false true], zbounds(2),...
                        node_coords(NOP_array(el,:),:))
                    bc = P_bcs{6};
                    el_nodes = NOP_array(el,:);

                    if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
                        b = one_el_neumann(6,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(6,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                end
                    % Now that we've gone through and applied all the
                    % Neumann conditions for this element, we check and see
                    % if there are any Dirichlet conditions that would
                    % override these Neumann conditions
                    
                    % Start with x boundaries
                    if onBoundary([true false false], xbounds(1),...
                        node_coords(NOP_array(el,:),:)) 
                        el_nodes = NOP_array(el,:);
                        % Is it a Drichlet condition? If so, then iterate by
                        % nodes. Otherwise, operate element-wise
                        bc = P_bcs{3};
                        if eq(bc{1},'D')
                            for node=el_nodes
                        
                                if boundaries(node, 3)
                                    % It's a Dirichlet B.C., so set the value
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue, 'function_handle') &&...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    elseif ~isa(bcvalue,'function_handle')
                                                b(node) = bcvalue(node);
                                                % If it's a function of
                                                % of the global node number
                                    end
                                
                                    items = find(Amatrix(node, :));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                            
                                
                                else
                           %         b(node) = 0;
                                    % This will override any Neumann
                                    % conditions implemented on nodes in
                                    % this element, even if they are not on
                                    % tihs boundary
                                end
                            end
                        end
                    elseif onBoundary([true false false],xbounds(2),...
                            node_coords(NOP_array(el,:),:))
                        el_nodes = NOP_array(el,:);
                        bc = P_bcs{4};
                        if eq(bc{1},'D')
                            for node=el_nodes
                                if boundaries(node,4)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle')&&...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    elseif ~isa(bcvalue, 'function_handle')
                                        b(node) = bcvalue(node);
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                        
                    end
                    
                    % Now check y-boundaries (columns one and 2)
                    if onBoundary([false true false],ybounds(1),...
                         node_coords(NOP_array(el,:),:))
                        bc = P_bcs{1};
                        el_nodes = NOP_array(el,:);
                        
                        if eq(bc{1},'D')
                            for node = el_nodes
                                if boundaries(node,1)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle')&&...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    else
                                        b(node) = bcvalue(node_coords(node,1));
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    elseif onBoundary([false true false],ybounds(2),...
                            node_coords(NOP_array(el,:),:))
                        bc = P_bcs{2};
                        el_nodes = NOP_array(el,:);
                        if eq(bc{1},'D')
                            for node = el_nodes
                                if boundaries(node,2)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle')&&...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    else
                                        b(node) = bcvalue(node_coords(node,1));
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    end
                    
                    % Last, z boundaries are checked for Dirichlet
                    % conditions
                    if onBoundary([false false true],zbounds(1),...
                            node_coords(NOP_array(el,:),:))
                        bc = P_bcs{5};
                        el_nodes = NOP_array(el,:);
                        if eq(bc{1},'D')
                            for node = el_nodes
                                if boundaries(node,5)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle')&&...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    else
                                        b(node) = bcvalue(node_coords(node,1));
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    elseif onBoundary([false false true],zbounds(2),...
                            node_coords(NOP_array(el,:),:))
                        bc = P_bcs{6};
                        el_nodes = NOP_array(el,:);
                        if eq(bc{1},'D')
                            bcvalue = bc{2};
                            for node = el_nodes
                                if boundaries(node,6)
                                    if ~isa(bcvalue,'function_handle')&&...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    else
                                        b(node) = bcvalue(node_coords(node,1));
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                                    
                            end
                        end
                    end
                
                el = el + 1;
            end
        end
    end
    
        
    
    bc2 = P_bcs{2};
    if eq(bc2{1}, 'N') &&  isa(bc2{2}, 'function_handle')
        [Amatrix, b] = Neumann_bc(Amatrix,b,Cvector,num_nodesx, num_nodesy,...
            num_nodesz, mu, kbm,bm_width,bc2, ybounds, xbounds, ...
            node_coords, issinglecap);
    end
    
    
 Pvector = Amatrix\b; % This returns P as a vector, ordered by nodes
error = 1;
limit = 0.001;
loop_counter = 1;

while error > limit
 
    % Now, track the concentration of two diffusing species
    %Cmatrix2 = spalloc(NNODES, NNODES, 6*NNODES);
    % The calculation of the diffusion term is the same as the loop from
    % Darcy's law -- included below, and with the addition of a term that comes
    % from convection
   
    Cmatrix1 = spalloc(NNODES, NNODES, 6*NNODES);
    Amatrix = spalloc(NNODES, NNODES, 6*NNODES);
    Cmatrix2 = spalloc(NNODES, NNODES, 6*NNODES);
    Cmatrix1 = calculate_concentration(Pvector,Cmatrix1,num_nodesx,...
        num_nodesy,num_nodesz, node_coords,kappa,mu,  Dalb);
    Cmatrix2 = calculate_concentration(Pvector,Cmatrix2,num_nodesx,num_nodesy,...
        num_nodesz,node_coords,kappa,mu, DIgA);
    % Create the boundary conditions vector
    b_concentration = zeros(NNODES,1);
    b_IgA = zeros(NNODES,1);
    
    for node=1:NNODES % Loop through all the nodes
        
        % Check to see if they are on any boundary
        if boundaries(node, 3)
            % For that boundary, determine what kind of bc we have, and
            % value of bc
            bc = C_bcs{3};
            bc2 = Ig_bcs{3};
            if (eq(bc{1},'D'))
                bcvalue = bc{2};
                if ~isa(bcvalue, 'function_handle')&&...
                                            isequal(size(bcvalue),[1,1])
                    b_concentration(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_concentration(node) = bcvalue(node);
                    
                end
                items = find(Cmatrix1(node, :));
                Cmatrix1(node,items) = 0;
                Cmatrix1(node,node) = 1.0;
            else
                if ~eq(bc{2}, 0)
                    %neumann_bc();
                end
            end
            if (eq(bc2{1}, 'D'))
                bcvalue = bc2{2};
                if ~isa(bcvalue, 'function_handle') && ...
                        isequal(size(bcvalue),[1,1])
                    b_IgA(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_IgA(node) = bcvalue(node);
                end
                items = find(Cmatrix2(node, :));
                Cmatrix2(node,items) = 0;
                Cmatrix2(node,node) = 1.0;
            end
        elseif boundaries(node,4)
            bc = C_bcs{4};
            bc2 = Ig_bcs{4};
            if eq(bc{1}, 'D')
                bcvalue = bc{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_concentration(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_concentration(node) = bcvalue(node);
                end
                items = find(Cmatrix1(node, :));
                Cmatrix1(node,items) = 0;
                Cmatrix1(node,node) = 1.0;
            end
            if (eq(bc2{1}, 'D'))
                bcvalue = bc2{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_IgA(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_IgA(node) = bcvalue(node);
                end
                items = find(Cmatrix2(node, :));
                Cmatrix2(node,items) = 0;
                Cmatrix2(node,node) = 1.0;
            end
            
        end
        if boundaries(node,1)
            bc = C_bcs{1};
            bc2 = Ig_bcs{1};
            if eq(bc{1},'D')
                bcvalue = bc{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_concentration(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_concentration(node) = bcvalue(node);
                end
                items = find(Cmatrix1(node,:));
            Cmatrix1(node,items) = 0;
            Cmatrix1(node,node) = 1;
            end
            if (eq(bc2{1}, 'D'))
                bcvalue = bc2{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_IgA(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_IgA(node) = bcvalue(node);
                end
                items = find(Cmatrix2(node, :));
                Cmatrix2(node,items) = 0;
                Cmatrix2(node,node) = 1.0;
            end
            
        elseif boundaries(node,2)
            bc = C_bcs{2};
            bc2 = Ig_bcs{2};
            if eq(bc{1},'D')
                bcvalue = bc{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_concentration(node) = bcvalue;
                end
                items = find(Cmatrix1(node,:));
            Cmatrix1(node,items) = 0;
            Cmatrix1(node,node) = 1;
            end
            if (eq(bc2{1}, 'D'))
                bcvalue = bc2{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_IgA(node) = bcvalue;
                end
                items = find(Cmatrix2(node, :));
                Cmatrix2(node,items) = 0;
                Cmatrix2(node,node) = 1.0;
            end
            
        end
        if boundaries(node,5)
            bc = C_bcs{5};
            bc2 = Ig_bcs{5};
            if eq(bc{1},'D')
                bcvalue = bc{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_concentration(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_concentration(node) = bcvalue(node);
                end
                items = find(Cmatrix1(node,:));
            Cmatrix1(node,items) = 0;
            Cmatrix1(node,node) = 1;
            end
            if (eq(bc2{1}, 'D'))
                bcvalue = bc2{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_IgA(node) = bcvalue;
                elseif ~isa(bcvalue, 'function_handle')
                    b_IgA(node) = bcvalue(node);
                end
                items = find(Cmatrix2(node, :));
                Cmatrix2(node,items) = 0;
                Cmatrix2(node,node) = 1.0;
            end
        elseif boundaries(node,6)
            bc = C_bcs{6};
            bc2 = Ig_bcs{6};
            if eq(bc{1},'D')
                bcvalue = bc{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_concentration(node) = bcvalue;
                else
%                     b(node) = bcvalue(Cmatrix1(node));
                end
                items = find(Cmatrix1(node,:));
            Cmatrix1(node,items) = 0;
            Cmatrix1(node,node) = 1;
            end
            if (eq(bc2{1}, 'D'))
                bcvalue = bc2{2};
                if ~isa(bcvalue, 'function_handle')&& ...
                        isequal(size(bcvalue),[1,1])
                    b_IgA(node) = bcvalue;
                end
                items = find(Cmatrix2(node, :));
                Cmatrix2(node,items) = 0;
                Cmatrix2(node,node) = 1.0;
            end
        end
        
        
    end
Cvector = Cmatrix1\b_concentration;
CIgvector = Cmatrix2\b_IgA;

Amatrix = calculate_pressure(Amatrix, num_nodesx,num_nodesy,num_nodesz,...
        node_coords,kappa,mu);
    b = zeros(NNODES,1);
    
     el = 1;
    for elx = 1:num_nodesx-1
        for ely = 1:num_nodesy-1
            for elz = 1:num_nodesz-1
                
                % First check the low x boundary (usually x = 0)
                if onBoundary([true false false], xbounds(1),...
                        node_coords(NOP_array(el,:),:)) 
                    el_nodes = NOP_array(el,:);
                    % Is it a Drichlet condition? If so, then iterate by
                    % nodes. Otherwise, operate element-wise
                    bc = P_bcs{3};
                    
                    if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
%                     
                        b = one_el_neumann(3,xbounds,boundaries,node_coords,...
                                    bc{2},b,num_nodesy,num_nodesx,...
                                    num_nodesz,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        
                         b = one_el_neumann(3,xbounds,boundaries,node_coords,...
                                    bc{2},b,num_nodesy,num_nodesx,...
                                    num_nodesz,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);

                    end
                    
                elseif onBoundary([true false false], xbounds(2),...
                        node_coords(NOP_array(el,:),:)) 
                    % That just checked boundary 4
                    % Is there a Dirichlet condition on this boundary?
                    bc = P_bcs{4};
                    el_nodes = NOP_array(el,:);
                    if eq(bc{1},'N') && (~isa(bc{2}, 'function_handle') && ~eq(bc{2},0))
                    
                        b = one_el_neumann(4,xbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(4,xbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                end
                
                % Now, check the y boundaries (1 & 2)
                 if onBoundary([false true false],ybounds(1),...
                         node_coords(NOP_array(el,:),:))
                     bc = P_bcs{1};

                        
                    if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
                        b = one_el_neumann(1,ybounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(1,ybounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                elseif onBoundary([false true false],ybounds(2),...
                        node_coords(NOP_array(el,:),:))
                     bc = P_bcs{2};
                     el_nodes = NOP_array(el,:);

                     if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0)
                          b = one_el_neumann(2,ybounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
%                      elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
%                          b = one_el_neumann(2,ybounds,boundaries,node_coords,...
%                             bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
%                             elx,ely,elz,Cvector,Amatrix...
%                                     ,kbm,bm_width,mu);
                     end
                end
                
                if onBoundary([false false true], zbounds(1),...
                        node_coords(NOP_array(el,:),:))
                    bc = P_bcs{5};
                    el_nodes = NOP_array(el,:);

                    if  eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0) 
                        b = one_el_neumann(5,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz...
                            ,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(5,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz...
                            ,elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                elseif onBoundary([false false true], zbounds(2),...
                        node_coords(NOP_array(el,:),:))
                    bc = P_bcs{6};
                    el_nodes = NOP_array(el,:);

                    if eq(bc{1},'N') && ~isa(bc{2}, 'function_handle') && ~eq(bc{2},0) 
                        b = one_el_neumann(6,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    elseif eq(bc{1},'N') && isa(bc{2},'function_handle')
                        b = one_el_neumann(6,zbounds,boundaries,node_coords,...
                            bc{2},b,num_nodesy,num_nodesx,num_nodesz,...
                            elx,ely,elz,Cvector,Amatrix...
                                    ,kbm,bm_width,mu);
                    end
                end
                    % Now that we've gone through and applied all the
                    % Neumann conditions for this element, we check and see
                    % if there are any Dirichlet conditions that would
                    % override these Neumann conditions
                    
                    % Start with x boundaries
                    if onBoundary([true false false], xbounds(1),...
                        node_coords(NOP_array(el,:),:)) 
                        el_nodes = NOP_array(el,:);
                        % Is it a Drichlet condition? If so, then iterate by
                        % nodes. Otherwise, operate element-wise
                        bc = P_bcs{3};
                        if eq(bc{1},'D')
                            for node=el_nodes
                        
                                if boundaries(node, 3)
                                    % It's a Dirichlet B.C., so set the value
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue, 'function_handle')&& ...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    elseif ~isa(bcvalue, 'function_handle')
                                        b(node) = bcvalue(node);
                                    end
                                
                                    items = find(Amatrix(node, :));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                            
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    elseif onBoundary([true false false],xbounds(2),...
                            node_coords(NOP_array(el,:),:))
                        el_nodes = NOP_array(el,:);
                        bc = P_bcs{4};
                        if eq(bc{1},'D')
                            for node=el_nodes
                                if boundaries(node,4)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle') && ...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    elseif ~isa(bcvalue, 'function_handle')
                                        b(node) = bcvalue(node);
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                        
                    end
                    
                    % Now check y-boundaries (columns one and 2)
                    if onBoundary([false true false],ybounds(1),...
                         node_coords(NOP_array(el,:),:))
                        bc = P_bcs{1};
                        el_nodes = NOP_array(el,:);
                        
                        if eq(bc{1},'D')
                            for node = el_nodes
                                if boundaries(node,1)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle') && ...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    elseif onBoundary([false true false],ybounds(2),...
                            node_coords(NOP_array(el,:),:))
                        bc = P_bcs{2};
                        el_nodes = NOP_array(el,:);
                        if eq(bc{1},'D')
                            for node = el_nodes
                                if boundaries(node,2)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle') && ...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    end
                    
                    % Last, z boundaries are checked for Dirichlet
                    % conditions
                    if onBoundary([false false true],zbounds(1),...
                            node_coords(NOP_array(el,:),:))
                        bc = P_bcs{5};
                        el_nodes = NOP_array(el,:);
                        if eq(bc{1},'D')
                            for node = el_nodes
                                if boundaries(node,5)
                                    bcvalue = bc{2};
                                    if ~isa(bcvalue,'function_handle') && ...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    elseif onBoundary([false false true],zbounds(2),...
                            node_coords(NOP_array(el,:),:))
                        bc = P_bcs{6};
                        el_nodes = NOP_array(el,:);
                        if eq(bc{1},'D')
                            bcvalue = bc{2};
                            for node = el_nodes
                                if boundaries(node,6)
                                    if ~isa(bcvalue,'function_handle') && ...
                                            isequal(size(bcvalue),[1,1])
                                        b(node) = bcvalue;
                                    end
                                    items = find(Amatrix(node,:));
                                    Amatrix(node,items) = 0;
                                    Amatrix(node,node) = 1.0;
                                else
%                                     b(node) = 0;
                                end
                            end
                        end
                    end
                
                el = el + 1;
            end
        end
    end
            
        
    bc2 = P_bcs{2};
    if eq(bc2{1}, 'N') &&  isa(bc2{2}, 'function_handle')
        %disp('Satisfied condition')
        [Amatrix, b] = Neumann_bc(Amatrix,b,Cvector,num_nodesx, ...
            num_nodesy, num_nodesz, mu, kbm,bm_width,bc2, ybounds, ...
            xbounds, node_coords, issinglecap);
    end
    Pvector2 = Amatrix\b;
    
    error = (Pvector2 - Pvector).^2;
    error = max(error);
    display(error);
    Pvector = Pvector2;
    
    loop_counter = loop_counter + 1;
end

Pmatrix = reshape(Pvector, num_nodesy*num_nodesz, num_nodesx);
Cmatrix = reshape(Cvector, num_nodesy*num_nodesz, num_nodesx);
CIgA = reshape(CIgvector, num_nodesy*num_nodesz, num_nodesx);
% Split for remaking as a 3d matrix
P = zeros(size(meshy));
C = zeros(size(meshy));
CIg = zeros(size(meshy));
count = 1;
for i=1:num_nodesz
    P(:,:,i) = Pmatrix(count:count+num_nodesy-1, :); 
   C(:,:,i) = Cmatrix(count:count+num_nodesy-1,:);
   CIg(:,:,i) = CIgA(count:count+num_nodesy-1,:);
    
    count = count+num_nodesy;
end
nout = max(nargout, 1) - 3;
if nout >= 1
    varargout{1} = Pvector;
end
if nout >= 2
    varargout{2} = Cvector;
end
if nout >= 3
   varargout{3} = b;
end
toc;
end

%  The order of the nodes is lower back
% left, upper back left, lower front left, upper front left, lower back
% right, upper back right, lower front right, upper front right. That is
% nodes are numbered first in y, then in z, then in x, for a right handed
% coordinate system

% This function returns the global node numbers of the 8 nodes in the given
% element elx, ely, elz

function [global_node_num] = element_nodes(ex, ey, ez, NNY, NNZ)

global_node_num = [0 0 0 0 0 0 0 0];
global_node_num(1) = (ex-1)*NNY*NNZ + (ez-1)*NNY + ey;
global_node_num(2) = (ex-1)*NNY*NNZ + (ez-1)*NNY + ey + 1;
global_node_num(3) = (ex-1)*NNY*NNZ + (ez-1)*NNY +ey + NNY;
global_node_num(4) = (ex-1)*NNY*NNZ + (ez-1)*NNY +ey + NNY + 1;
global_node_num(5) = (ex-1)*NNY*NNZ + (ez-1)*NNY + NNZ*NNY + ey;
global_node_num(6) = global_node_num(5) + 1;
global_node_num(7) = global_node_num(5) + NNY;
global_node_num(8) = global_node_num(6) + NNY;
end

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
% are xi, etta, zeta coordinate and weight. These map to an element that goes
% from -1 to 1

% This needs to change to be a function of the boundary
gauss_points2 = zeros(4,4);
if boundary == 1 
    gauss_points2(1,:) = [-0.5774 -1 -0.5774 1];
    gauss_points2(2,:) = [-0.5774 -1 -0.5774 1];
    gauss_points2(3,:) = [-0.5774 -1 0.5774 1];
    gauss_points2(4,:) = [-0.5774 -1 0.5774 1];
    gauss_points2(5,:) = [0.5774 -1 -0.5774 1];
    gauss_points2(6,:) = [0.5774 -1 -0.5774 1];
    gauss_points2(7,:) = [0.5774 -1 0.5774 1];
    gauss_points2(8,:) = [0.5774 -1 0.5774 1];
elseif boundary == 2
   gauss_points2(1,:) = [-0.5774 1 -0.5774 1];
    gauss_points2(2,:) = [-0.5774 1 -0.5774 1];
    gauss_points2(3,:) = [-0.5774 1 0.5774 1];
    gauss_points2(4,:) = [-0.5774 1 0.5774 1];
    gauss_points2(5,:) = [0.5774 1 -0.5774 1];
    gauss_points2(6,:) = [0.5774 1 -0.5774 1];
    gauss_points2(7,:) = [0.5774 1 0.5774 1];
    gauss_points2(8,:) = [0.5774 1 0.5774 1];
elseif boundary == 3
    gauss_points2(1,:) = [-1 -0.5774 -0.5774 1];
    gauss_points2(2,:) = [-1 0.5774 -0.5774 1];
    gauss_points2(3,:) = [-1 -0.5774 0.5774 1];
    gauss_points2(4,:) = [-1 0.5774 0.5774 1];
    gauss_points2(5,:) = [-1 -0.5774 -0.5774 1];
    gauss_points2(6,:) = [-1 0.5774 -0.5774 1];
    gauss_points2(7,:) = [-1 -0.5774 0.5774 1];
    gauss_points2(8,:) = [-1 0.5774 0.5774 1];
elseif boundary == 4
    gauss_points2(1,:) = [1 -0.5774 -0.5774 1];
    gauss_points2(2,:) = [1 0.5774 -0.5774 1];
    gauss_points2(3,:) = [1 -0.5774 0.5774 1];
    gauss_points2(4,:) = [1 0.5774 0.5774 1];
    gauss_points2(5,:) = [1 -0.5774 -0.5774 1];
    gauss_points2(6,:) = [1 0.5774 -0.5774 1];
    gauss_points2(7,:) = [1 -0.5774 0.5774 1];
    gauss_points2(8,:) = [1 0.5774 0.5774 1];
elseif boundary == 5
    gauss_points2(1,:) = [-0.5774 -0.5774 -1 1];
    gauss_points2(2,:) = [-0.5774 0.5774 -1 1];
    gauss_points2(3,:) = [-0.5774 -0.5774 -1 1];
    gauss_points2(4,:) = [-0.5774 0.5774 -1 1];
    gauss_points2(5,:) = [0.5774 -0.5774 -1 1];
    gauss_points2(6,:) = [0.5774 0.5774 -1  1];
    gauss_points2(7,:) = [0.5774 -0.5774 -1 1];
    gauss_points2(8,:) = [0.5774 0.5774 -1  1];
else
    gauss_points2(1,:) = [-0.5774 -0.5774 1 1];
    gauss_points2(2,:) = [-0.5774 0.5774 1 1];
    gauss_points2(3,:) = [-0.5774 -0.5774 1 1];
    gauss_points2(4,:) = [-0.5774 0.5774 1 1];
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
% Implement a Neumann boundary by adding the boundary integral to the
% b vector (passed in to this function as bvalues

% This applies the Neumann boundary condition to the entire boundary. USing
% the way the convection-diffusion code is structured now, we just want to
% do that on the particular element we have now

for ely=1:nny-1
    for elx=1:nnx-1
        for elz=1:nnz-1
            % get node coords for this elemnet
            nodes = element_nodes(elx,ely,elz, nny, nnz);
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
        end
    end
end

b = bvalues;
end

% This function returns a 1 (true) if the given element or node is on the given
% boundary, and 0 (false) otherwise. node_or_element_coords contains either
% just three values (x,y,z coordinates of one node), or it can contain the  
% x, y, z coordinates of 8 nodes, if we are checking elements. Direction is
% a logical 3 vector of the form [x y z]. It has a 1 in the position
% showing which direction boundary we are testing. Value is then the value
% for that coordinate at this boundary

function [yes_no] = onBoundary(direction, value, node_or_element_coords)
[rows, cols] = size(node_or_element_coords);

if rows == 1
    yes_no = (node_or_element_coords(direction) == value);
else
    count = 1;
    while count <= rows && node_or_element_coords(count,direction) ~= value
        count = count + 1;
    end
    if count ~= 9
       yes_no = (node_or_element_coords(count,direction) == value);
    else
        yes_no = false;
    end
end

end
