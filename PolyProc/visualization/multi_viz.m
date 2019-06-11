function [adj_p_coord] = multi_viz(gid_map, modal_2, varargin)
% grain_pmdf vizualizes grain boundary colored by the number of adjacent particles
%==========================================================================
% FILENAME:          grain_pmdf.m
% DATE:              1 May, 2019        
% PURPOSE:           visualize grain boundaries together with particles
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%    modal_2   : (array) n*3 array containing x,y,z coordinates
%OPTIONAL :
%    goi       : (1*n array) list of grain ID to investigate
%                   default = entire list of grains
%    distance_threshold : (double) distance threshold to count number of particles
%                         default = 2 voxels
%
%OUT :
%    3D visualization of grain boundadry colored by adjacent particles density
%    adj_p_coord : (array) m*3 array containing x,y,z coordinates of adj features of modal_2 
%==========================================================================
%EXAMPLE :
%    - investigate entire grains
%       multi_viz(gid_map_1, P_position_t0)
%
%    - investigate specified grains with default distance threshold '2'
%       multi_viz(gid_map_1, P_position_t0, 'goi', [5 10 12 20])
%
%    - investigate specified grains with user-specified distance threshold '1'
%       multi_viz(gid_map_1, P_position_t0, 'goi', [5 10 12 20], 'distance_threshold', 1)

%==========================================================================

if any(strcmp(varargin,'goi')) 
    idx = find(strcmp(varargin,'goi'))+1;
    goi = varargin{idx};
    
elseif ~any(strcmp(varargin,'goi'))
    goi = unique(gid_map);
    goi = goi(2:end);
end

if any(strcmp(varargin,'distance_threshold')) 
    distance_threshold = find(strcmp(varargin,'distance_threshold'))+1;
else 
    distance_threshold = 2; % default distance_threshold
end

figure
fprintf('Searching for nearby features and plotting them ...\n');

    for i = 1:length(goi)
        FV = isosurface(gid_map == goi(i));
        
        if size(FV.vertices,1) > 10
            patch(smoothpatch(FV,1,3), ...
                        'facecolor', 'b', ...
                        'edgecolor', 'none', ...
                        'facealpha', 0.3);
            hold on


            [~,D] = knnsearch(FV.vertices,modal_2);
            adj_p = D<distance_threshold;

            scatter3(modal_2(adj_p,1),modal_2(adj_p,2),modal_2(adj_p,3),12,'filled','MarkerFaceColor','r');

        end
        
        hold on
    end

view([30,20]); axis equal; camlight left; camlight right; material metal;
axis equal
axis off
legend('Grain','Secondary feature')
adj_p_coord = modal_2(adj_p,:);
hold on

end
