function [ gauss_points ] = get_gauss_points( )
%#codgen
%get_gauss_points Returns a matrix of the gauss points in order, with col 1
%as the x coordinate, col 2 y coordinate, col 3 z coordinate, and col 4
%gauss point weight
%   Detailed explanation goes here

num_gp = 8;
gp = 1/sqrt(3);
gauss_points = zeros(num_gp,4);
    gauss_points(1,:) = [-gp -gp, -gp  1];
    gauss_points(2,:) = [-gp gp -gp 1];
    gauss_points(3,:) = [-gp -gp gp 1];
    gauss_points(4,:) = [-gp gp  gp 1];
    gauss_points(5,:) = [gp -gp, -gp  1];
    gauss_points(6,:) = [gp gp, -gp  1];
    gauss_points(7,:) = [gp -gp, gp  1];
    gauss_points(8,:) = [gp gp, gp  1];
end

