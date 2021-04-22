% MATLAB R2018a
%
% FUNCTION NAME:
%   read_vtk
%
% DESCRIPTION:
%   Read a .vtk mesh file for plotting
%
% INPUT:
%   filename - (string) loaction of the file to load
%
% OUTPUT:
%   surf - (struct) a structure with triangle IDs and vertex coordinates
%   header - (struct) a structure with misc information about the .vtk file
%   Side effects: none
%
% ASSUMPTIONS AND LIMITATIONS:
%   The .vtk file must be of a closed triangle mesh
%
function [surf, header] = read_vtk(filename)

% attempt to open the vtk file
fid = fopen(filename, 'rb');

if (fid < 0)
    fprintf('could not open file %s\n', filename);
    return
end

% read the .vtk header information
str = fgetl(fid);
header.Filename = filename;
header.Format = str(3:5);
header.Version = str(end-2:end);
header.Header = fgetl(fid);
header.DatasetFormat = lower(fgetl(fid));

% read the data
str = lower(fgetl(fid));
header.DatasetType = str(9:end);

if ~strcmpi(header.DatasetType, 'polydata')
    fprintf('ERROR: only the "polydata" vtk data type is supported');
    return
end

while ~feof(fid)
    str = lower(fgetl(fid));
    
    % read point data (coordinates)
    if (contains(str, 'points'))   
        % find the number of vertices in the file
        tmp = textscan(str,'%s %d %s');
        npoints = tmp{2};
        
        surf.npoints = npoints;    
        [values,count] = fscanf(fid,'%f %f %f',3*npoints);

        if count ~= 3*npoints
            fprintf('ERROR: there was a problem reading the vertices');
            return
        end
        
        values = reshape(values, 3, count/3);
        surf.vtx = values;
        
    % read polygon data (triangles)
    elseif (contains(str, 'polygons'))
        % find the number of vertices in the file
        tmp = textscan(str,'%s %d %d');
        ntriangles = tmp{3};
        
        surf.ntriangles = ntriangles / 4;
        [values, count] = fscanf(fid, '%d %d %d %d', 4*ntriangles);
       
        if count ~= ntriangles
            fprintf('ERROR: there was a problem reading the triangles');
            return
        end
        
        values = reshape(values, 4, count/4);
        surf.tri = values(2:4,:)' + 1;
    end
end

end