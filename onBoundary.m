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
