function [serialized_map, id_mapping] = serialize_grain_ids(grain_map)
% serialize_grain_ids - Serializes the grain IDs in a 3D grain map.
%
% PURPOSE: Convert the grain IDs in a 3D map to a continuous range of IDs,
% preserving zero values, which represent background or no-grain areas.
%
% INPUT:
%   grain_map - 3D numeric array where each unique non-zero integer 
%               represents a distinct grain ID.
%
% OUTPUT:
%   serialized_map - 3D array of the same size as grain_map with serialized
%                    grain IDs, where IDs are continuous starting from 1.
%   id_mapping     - Struct mapping original grain IDs to serialized IDs.
%
% EXAMPLE USAGE:
%   [serialized_map, id_mapping] = serialize_grain_ids(grain_map);

% Initialize output map
serialized_map = zeros(size(grain_map), 'int32');

% Find unique, non-zero grain IDs in the map
unique_ids = unique(grain_map);
unique_ids(unique_ids == 0) = [];  % Exclude zero

% Create a mapping of old IDs to new serialized IDs
new_ids = 1:length(unique_ids);
id_mapping = containers.Map(unique_ids, new_ids);

% Serialize grain IDs in the map
for i = 1:length(unique_ids)
    serialized_map(grain_map == unique_ids(i)) = new_ids(i);
end

% Return id_mapping as a struct for easy access
id_mapping = struct('original_ids', unique_ids, 'serialized_ids', new_ids);

end
