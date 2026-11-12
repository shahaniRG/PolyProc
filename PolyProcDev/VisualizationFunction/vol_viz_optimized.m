function [] = vol_viz_optimized(gid_map, numElement, varargin)
% Highly Optimized GPU-Accelerated 3D Volume Visualization

% Move data to GPU for speedup
if ~isa(gid_map, 'gpuArray')
    gid_map = gpuArray(gid_map);
end

% Extract unique grain IDs
gid_list = unique(gid_map(:));
gid_list = gid_list(gid_list > 0);  % Exclude background (0)
num_grains = length(gid_list);

% Initialize color list
color_gid_list = rand(num_grains, 3, 'gpuArray'); % Default random coloring

% Handle optional arguments
selective = [];
exclude = [];
if any(strcmp(varargin, 'selective'))
    idx = find(strcmp(varargin, 'selective')) + 1;
    selective = varargin{idx};
elseif any(strcmp(varargin, 'exclude'))
    idx = find(strcmp(varargin, 'exclude')) + 1;
    exclude = varargin{idx};
end

% Filter gid_list based on selection/exclusion
if ~isempty(selective)
    gid_list = intersect(gid_list, gpuArray(selective));
elseif ~isempty(exclude)
    gid_list = setdiff(gid_list, gpuArray(exclude));
end

% Efficient coloring definition
fprintf('Defining coloring on GPU...\n');
if any(strcmp(varargin, 'orientation'))
    % Orientation-based coloring
    idx = find(strcmp(varargin, 'orientation')) + 1;
    grain_rodV = gpuArray(varargin{idx});
    crystal = varargin{idx + 1};
    cs = crystalSymmetry(crystal);
    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;

    % Vectorized orientation coloring
    grain_orientation = arrayfun(@(i) rod2quat(grain_rodV(i, :)), 1:num_grains, 'UniformOutput', false);
    color_gid_list = cell2mat(arrayfun(@(i) ...
        oM.orientation2color(orientation('quaternion',grain_orientation{i}, cs)), ...
        1:num_grains, 'UniformOutput', false));

elseif any(strcmp(varargin, 'completness'))
    % Completeness-based coloring
    idx = find(strcmp(varargin, 'completness')) + 1;
    comp = gpuArray(varargin{idx});
    comp_mean = arrayfun(@(gid) mean(comp(gid_map == gid)), gid_list);
    comp_color = jet(101);
    color_gid_list = comp_color(ceil((comp_mean - min(comp_mean)) / ...
        (max(comp_mean) - min(comp_mean) + eps) * 100), :);

elseif any(strcmp(varargin, 'size'))
    % Size-based coloring
    idx = find(strcmp(varargin, 'size')) + 1;
    numElement = gpuArray(varargin{idx});
    size_color = jet(max(numElement(:, 2)));
    color_gid_list = size_color(numElement(:, 2), :);

elseif any(strcmp(varargin, 'neighbor'))
    % Neighbor-based coloring
    adj = gpuArray(imRAG(gid_map));
    adj = adj(adj(:, 1) < adj(:, 2), :);  % Avoid duplicate pairs
    [~, F] = mode(adj(:));
    neighbor_color = jet(F);

    % Vectorized neighbor calculation
    numNeighbor = arrayfun(@(gid) sum(sum(adj == gid)), gid_list);
    color_gid_list = neighbor_color(max(1, min(numNeighbor, F)), :);
end

% Visualization
figure;
%parpool('local'); % Ensure parallel pool is active
for i = 1:num_grains
    % Create isosurface on GPU
    FV = isosurface(gather(gid_map == gid_list(i)));  % Gather only necessary data
    if ~isempty(FV.vertices)
        % Plotting with minimal data transfer
        patch('Vertices', FV.vertices, 'Faces', FV.faces, ...
              'FaceColor', gather(color_gid_list(i, :)), ...
              'EdgeColor', 'none', 'FaceAlpha', 1);
        hold on;
    end
end

% Plot settings
view([30, 20]);
axis equal;
material dull;
camlight headlight;

% Display colorbar based on coloring method
if any(strcmp(varargin, 'completness'))
    colormap(comp_color);
    caxis([min(comp_mean) max(comp_mean)]);
    colorbar;
elseif any(strcmp(varargin, 'size'))
    colormap(size_color);
    caxis([0 max(numElement(:, 2))]);
    colorbar;
elseif any(strcmp(varargin, 'neighbor'))
    colormap(neighbor_color);
    caxis([0 F]);
    colorbar;
end

end
