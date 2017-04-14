function [ meshgrid_matrix ] = make_matrix( node_vector, x,y,z )
% This function takes in a vector in node order, that is with values
% associated with nodes by row number, and reorders it to be in a meshgrid
% format, where the nodes are numbered first in y, then in z, and last in
% x. x,y, z are the meshgrid outputs

numx = size(x,2);
numy = size(y,1);
numz = size(z,3);

meshgrid_matrix = zeros(size(y));

node_matrix = reshape(node_vector,numy*numz,numx);
count = 1;
for i=1:numz
    meshgrid_matrix(:,:,i) = node_matrix(count:count + numy-1,:);
    count = count + numy;
end

end

