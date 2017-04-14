function [ node_coords ] = get_node_coords( meshx, meshy, meshz)
%get_node_coords This function returns a matrix, node_coords, that gives
%the real space x,y, and z coordinates of the node i, in row i of the
%matrix
%   This is a matrix with N rows, where N is the total number of nodes, and
%   3 columns. Column 1 is the x-coordinate, col 2 is the y-coordinate, and
%   col 3 is the z-coordinate of the point number i, where i is the row
%   number in node_coords

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

node_coords = zeros(num_nodesx*num_nodesy*num_nodesz, 3);
NNODES = num_nodesx*num_nodesy*num_nodesz;
for i=1:NNODES
    node_coords(i,1:3) = [meshx(i) meshy(i) meshz(i)];
end

end

