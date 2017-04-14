function [ new_vector ] = convert_data(old_vector,x,y,z )
% This function converts data stored in matrix with the order y, x, z to a
% matrix stored with the order y, z, x. x, y,z should be from meshgrid

num_ynodes = size(y,1);
num_xnodes = size(x,2);
num_znodes = size(z,3);

new_vector = zeros(size(old_vector));
zloops = 1;
for xcounter=1:num_xnodes
    for zcounter=1:num_znodes
        new_vector((xcounter-1)*num_ynodes*num_znodes+...
            (zcounter-1)*num_ynodes+1:...
            (zcounter)*num_ynodes+(xcounter-1)*num_ynodes*num_znodes,:) = ...
            old_vector((zcounter-1)*num_xnodes*...
            num_ynodes+(xcounter-1)*num_ynodes+1:...
            (zcounter-1)*num_xnodes*num_ynodes+xcounter*num_ynodes,:);
        zloops = zloops + 1;
    end
end

end

