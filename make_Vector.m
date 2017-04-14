% Convert a matrix in meshgrid form to the right vector form

function [vector] = make_Vector(meshgrid_m)

dims = size(meshgrid_m);
num_nodesy = dims(1);
num_nodesz = dims(3);
num_nodesx = dims(2);

array = zeros(num_nodesy, num_nodesz);
count = 0;

vector = zeros(num_nodesy*num_nodesz*num_nodesx,1);
for i=1:num_nodesx
    array(1:num_nodesy,1:num_nodesz) = meshgrid_m(:,i,:);
    vector(count+1:num_nodesy*num_nodesz*i) = ...
        reshape(array,num_nodesy*num_nodesz,1);
    count = i*num_nodesy*num_nodesz;
    
end


end