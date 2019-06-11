function [neighbor_angle, neighbor_angle_axis, neighbor_surface] = ...
                    boundary_viz(gid_map, grain_rodV, crystal, adj, varargin)
% boundary_visualization visualize neighboring grains of specified grain
%==========================================================================
% FILENAME:          boundary_visualization.m
% DATE:              1 May, 2019        
% PURPOSE:           visualize boundary between neighboring grains
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%    grain_rodV: n*3 array with rodrigues vectors of grains
%    crystal   : (string) 'cubic' for cubic crystal structure
%    adj       : n*2 array shows grain adjacency in number ascending order
%
%OPTIONAL :
%    goi       : (1*n array) list of grain ID to investigate
%                   default = entire list of grains
%
%OUT :
%    3D reconstruction of grain boundary with neighboring grains
%       ( grain color is based on normal direction of grain)
%    neighbor_angle     : n(number of neighboring grains)*1 array with misorientation 
%    neighbor_angle_axis: n(number of neighboring grains)*3 array with tilt axis    
%    neighbor_surface   : n(number of neighboring grains)*1 array with grain boundary area (in voxel^2)
%==========================================================================
%EXAMPLE :
%    [neighbor_angle_2, neighbor_angle_axis_2, neighbor_surface_2]= ...
%         boundary_viz(gid_map_2, grain_rodV_2, 'cubic', adj_2,'goi', 20);
%==========================================================================
cs = crystalSymmetry(crystal);

if any(strcmp(varargin,'goi')) 
    idx = find(strcmp(varargin,'goi'))+1;
    goi = varargin{idx};
    
elseif ~any(strcmp(varargin,'goi'))
    goi = unique(gid_map);
    goi = goi(2:end);
end

neighbor_surface = zeros(1,length(goi));
neighbor_angle = zeros(1,length(goi));
neighbor_angle_axis = cell(1,length(goi));

gid_list = unique(gid_map(:));
gid_list = gid_list(2:end);

figure; 
fprintf('Calculating misorientations and plotting grain boundaries ...\n');
for j = 1:length(goi)
    
            if any(grain_rodV(gid_list==goi(j),:))
                q1 = rodrigues2quat(vector3d(grain_rodV(gid_list==goi(j),:)));
            else
                sprintf('Grain %d does not exist in matrix\n',goi(j))
            end
        o1 = orientation(q1, cs);

        % Find its first-order nearest neighbors
        spec_gid_list = sum( adj( adj(:,1) == goi(j) | adj(:,2) == goi(j),:), 2) - double(goi(j));

        
        for i = 1:numel(spec_gid_list)

            q2 = rodrigues2quat(vector3d(grain_rodV(gid_list==spec_gid_list(i),:)));
            o2 = orientation(q2, cs);

            gidK = gid_map ~= goi(j) & gid_map ~= spec_gid_list(i);
            gidLocal = double(gid_map);
            gidLocal(gidK) = NaN;

            % Create the isosurface of gidLocal
            % Smooth the mesh
            FV = isosurface(gidLocal, (99*goi(j)+double(spec_gid_list(i)))/100);

            if size(FV.vertices,1) > 10

                % Plot the mesh of the grain boundary colored according to
                % misorientation angle between adjacent grains
                patch(smoothpatch(FV,1,3), ...
                    'facecolor', 'flat', ...
                    'edgecolor', 'none', ...
                    'facevertexcdata', double(angle(o1,o2)/degree), ...
                    'facealpha', 0.75);
                hold on;
            
                FV2 = smoothpatch(FV,1,3);
                verts = FV2.vertices;
                faces = FV2.faces;
                e21 =  verts(faces(:, 2), :) - verts(faces(:, 1), :);
                e31 = verts(faces(:, 3), :) - verts(faces(:, 1), :);
                each_area = cross(e21, e31, 2);
                neighbor_surface(i,j) = 1/2 * sum(sqrt(sum(each_area.^2, 2)));
                neighbor_angle(i,j) = double(angle(o1,o2)/degree);
                ax = axis(o1,o2);
                hkl = round(Miller(ax.x,ax.y,ax.z,cs),'maxHKL', 6);
                neighbor_angle_axis{i,j} = round([hkl.hkl]);
                
            else 
                warning('Boundary between grains %d and %d does not have enough vertices to display\n ', j, spec_gid_list(i));
            end
        end

end

view([30,20]); axis equal; axis off;
c = colorbar;
ylabel(c, 'Misorientation (deg)');
colormap(jet)

end
