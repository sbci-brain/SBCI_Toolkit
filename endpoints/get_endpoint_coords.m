% MATLAB R2018a
%
% FUNCTION NAME:
%   get_endpoint_coords
%
% DESCRIPTION:
%   Get the cartesian coordinates of the given endpoints on the given mesh 
%
% INPUT:
%   triangulation - (struct) The mesh to calculate over
%   intersections - (struct) The intersection information for a subject
%
% OUTPUT:
%   coords_in - (matrix) The cartesian coordinates for one end of each fiber
%   coords_out - (matrix) The cartesian coordinates for one end of each fiber
%   surf_in - (matrix) Returns a 0 for left hemisphere and 1 for right
%   hemisphere for each fiber
%   surf_out - (matrix) Returns a 0 for left hemisphere and 1 for right
%   hemisphere for each fiber
%
% ASSUMPTIONS AND LIMITATIONS:
%   The intersection and triangulation must come from the same resolution. 
%

function [coords_in, coords_out, surf_in, surf_out] = get_endpoint_coords(triangulation, intersections)

lh_V = triangulation.lh_V;
rh_V = triangulation.rh_V;
lh_T = triangulation.lh_T;
rh_T = triangulation.rh_T;

surf_in = intersections.surf_in;
surf_out = intersections.surf_out;
tri_in = intersections.tri_in;
tri_out = intersections.tri_out;
pt_in = intersections.pt_in;
pt_out = intersections.pt_out;

% Convert points from barycentric to cartesian
lh_in_points = pt_in(surf_in == 0,:);
lh_in_vertices = reshape(lh_V(lh_T(tri_in(surf_in == 0),:), :), [], 3, 3); 
lh_in_cartesian = squeeze(sum(lh_in_vertices .* lh_in_points, 2));

% Convert points from barycentric to cartesian
lh_out_points = pt_out(surf_out == 0,:);
lh_out_vertices = reshape(lh_V(lh_T(tri_out(surf_out == 0),:), :), [], 3, 3); 
lh_out_cartesian = squeeze(sum(lh_out_vertices .* lh_out_points, 2));
  
% Convert points from barycentric to cartesian
rh_in_points = pt_in(surf_in == 1,:);
rh_in_vertices = reshape(rh_V(rh_T(tri_in(surf_in == 1),:), :), [], 3, 3); 
rh_in_cartesian = squeeze(sum(rh_in_vertices .* rh_in_points, 2));

% Convert points from barycentric to cartesian
rh_out_points = pt_out(surf_out == 1,:);
rh_out_vertices = reshape(rh_V(rh_T(tri_out(surf_out == 1),:), :), [], 3, 3); 
rh_out_cartesian = squeeze(sum(rh_out_vertices .* rh_out_points, 2));

% Concatenate all the coordinates
coords_in = zeros(size(pt_in));
coords_in(surf_in == 0,:) = lh_in_cartesian;
coords_in(surf_in == 1,:) = rh_in_cartesian;

% Concatenate all the coordinates
coords_out = zeros(size(pt_out));
coords_out(surf_out == 0,:) = lh_out_cartesian;
coords_out(surf_out == 1,:) = rh_out_cartesian;

end