function [neighbor_angle, neighbor_angle_axis, neighbor_surface, p_storage] = ...
                    pmdf(gid_map, grain_rodV, crystal, adj, modal_2, varargin)
% pmdf_single
%==========================================================================
% FILENAME:          boundary_visualization.m
% DATE:              1 May, 2019        
% PURPOSE:           visualize neighboring grains
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%    grain_rodV: n*3 array with rodrigues vectors of grains
%    crystal   : (string) 'cubic' for cubic crystal structure
%    adj       : n*2 array shows grain adjacency in number ascending order
%
%OPTIONAL :
%    goi       : grain of interest
%                  default = entire grains
%    distance_threshold : (double) distance threshold to count number of particles
%                         default = 2 voxels
%
%OUT :
%    neighbor_angle     : n(number of neighboring grains)*1 array with misorientation 
%    neighbor_angle_axis: n(number of neighboring grains)*3 array with tilt axis    
%    neighbor_surface   : n(number of neighboring grains)*1 array with grain boundary area (in voxel^2)
%    p_storage          : n(number of neighboring grains)*1 array with the number of adjacent particles 
%==========================================================================
%EXAMPLE :
%    - PMDF for entire grains
%    [neighbor_angle_1, neighbor_angle_axis_1, neighbor_surface_1, p_storage]= ...
%               pmdf(gid_map_1, grain_rodV_1, 'cubic', adj_1, P_position_t0);
%
%    - PMDF for single grain
%    [neighbor_angle_1, neighbor_angle_axis_1, neighbor_surface_1, p_storage]= ...
%               pmdf(gid_map_1, grain_rodV_1, 'cubic', adj_1, P_position_t0, 'goi', 117);
%==========================================================================

if any(strcmp(varargin,'goi')) 
    idx = find(strcmp(varargin,'goi'))+1;
    goi = varargin{idx};
    
elseif ~any(strcmp(varargin,'goi'))
    goi = unique(gid_map);
    goi = goi(2:end);
end



if any(strcmp(varargin,'distance_thrshold')) 
    distance_threshold = find(strcmp(varargin,'distance_thrshold'))+1;
else 
    distance_threshold = 2; % default distance_threshold
end

n = 1;
p_storage = zeros(1,length(goi));

cs = crystalSymmetry(crystal);

neighbor_surface = zeros(1,length(goi));
neighbor_angle = zeros(1,length(goi));
neighbor_angle_axis = cell(1,length(goi));

redundancy = [];

for j = 1:length(goi)

    gid_list = unique(gid_map(:));
    gid_list = gid_list(2:end);
        if any(grain_rodV(gid_list==goi(j),:))
            q1 = rodrigues2quat(vector3d(grain_rodV(gid_list==goi(j),:)));
        else
            sprintf('Grain %d does not exist in matrix',goi(j))
        end
    o1 = orientation(q1, cs);

    % Find its first-order nearest neighbors
    spec_gid_list = sum( adj( adj(:,1) == goi(j) | adj(:,2) == goi(j),:), 2) - double(goi(j));

    for i = 1:numel(spec_gid_list)

        redundancy_check = sum((redundancy == goi(j)) + (redundancy == spec_gid_list(i)),2 );
        if ~any(ismember(redundancy_check,2))
            
            redundancy = [redundancy ; goi(j) spec_gid_list(i)];

            q2 = rodrigues2quat(vector3d(grain_rodV(gid_list==spec_gid_list(i),:)));
            o2 = orientation(q2, cs);

            gidK = gid_map ~= goi(j) & gid_map ~= spec_gid_list(i);
            gidLocal = double(gid_map);
            gidLocal(gidK) = NaN;

            FV = isosurface(gidLocal, (99*goi(j)+double(spec_gid_list(i)))/100);
            FV = smoothpatch(FV,1,3);

                if size(FV.vertices,1) > 10

                        [~,D] = knnsearch(FV.vertices,modal_2);
                        adj_p = D<distance_threshold;
                        p_storage(i,j) = sum(adj_p);
                        n = n+1;

                    verts = FV.vertices;
                    faces = FV.faces;
                    e21 =  verts(faces(:, 2), :) - verts(faces(:, 1), :);
                    e31 = verts(faces(:, 3), :) - verts(faces(:, 1), :);
                    each_area = cross(e21, e31, 2);
                    neighbor_surface(i,j) = 1/2 * sum(sqrt(sum(each_area.^2, 2)));
                    neighbor_angle(i,j) = double(angle(o1,o2)/degree);
                    ax = axis(o1,o2);
                    hkl = round(Miller(ax.x,ax.y,ax.z,cs),'maxHKL', 6);
                    neighbor_angle_axis{i,j} = round([hkl.hkl]);

                    hold on;
                else 
                    warning('Boundary between grains %d and %d does not have enough vertices to display\n ', goi(j) ,spec_gid_list(i));
                end
        else
            fprintf('pmdf analysis for grain boundary between %d and %d is alreay done.\n ', goi(j) ,spec_gid_list(i));
        end

    end
end

%% pmdf_figure


x = 0:2.5:60;
x_plot = 1.25:2.5:58.75;
y_particle = [];
y_surface = [];

for i = 1:length(x)-1
    
    idx = neighbor_angle>x(i) & neighbor_angle<x(i+1);
    
    y_surface(i) = sum(sum(sum (neighbor_surface .* idx)));
    y_particle(i) = sum(sum(sum  ( p_storage .* idx)));
    
end

y_particle_norm = y_particle ./ sum(sum(sum(y_particle)));
y_surface_norm = y_surface ./ sum(sum(sum(y_surface)));

yyaxis left
particle = bar(x_plot,y_particle_norm, 'FaceColor',[0, 0.4470, 0.7410]);
particle.BarWidth = 0.5;
ylabel('Vol. Fraction of Particles','fontsize',12,'fontname', 'Arial')

hold on
yyaxis right
bar(x_plot,y_surface_norm,'FaceColor', [0.8500, 0.3250, 0.0980], 'FaceAlpha',0.5);
ylabel('Area Fraction of Grain Boundaries','fontsize',12,'fontname', 'Arial')

s = ['Misorientation (' char(176) ')'];
xlabel(s,'fontsize',12,'fontname', 'Arial')

% led = legend;
% led.Location = 'northwest';
% led.String = {'Vol. Frac. Particles','Area Frac. G.B.'};
legend('Vol. Frac. Particles','Area Frac. G.B.','locattion', 'northwest');

end
