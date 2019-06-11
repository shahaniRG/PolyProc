function [] = vol_viz(gid_map, varargin)
% volume_visualization reconstruct 3d volume
%==========================================================================
% FILENAME:          volume_visualization.m
% DATE:              1 May, 2019        
% PURPOSE:           reconstruct 3d volume
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%
%OPTIONAL 1 : specify grain coloring criterion 
%     |---no specification  - random coloring (default)
%
%     |---'orientation'     - (string) coloring based on crystallographic dirction of grain
%                             (grain_rodV & crystal structure should be input)
%     |---'size'            - (string) coloring based on grain volume
%                             (numElement should be input)
%     |---'completeness'    - (string) coloring based on completeness of grain
%                             (comp should be input)
%     |---'neighbor'        - (string) coloring based on the number of neighboring grain
%                             
%OPTIONAL 2 : entire grain vs specific griain
%     |---no specification  - entire grains (default)
%
%     |---'selective'   - (1*n array) list of grains id to visualize only
%
%     |---'exclude'     - (1*n array) list of grains id to exclude
%
%OUT :
%    3D reconstructed object
%==========================================================================
% example
% 0. no specification
% vol_viz(gid_map_1)
%
% 1. orientation
% vol_viz(gid_map_1, 'orientation', grain_rodV_1, 'cubic')
% vol_viz(gid_map_2, 'orientation',grain_rodV_2, 'cubic', 'selective', [1 5 3 7])
% 
% 2. size
% vol_viz(gid_map_1, 'size', numElement_1)
% vol_viz(gid_map_2, 'size', numElement_2, 'selective', [1 5 3 7])
% 
% 3. completeness
% vol_viz(gid_map_1, 'completness', comp_1)
% vol_viz(gid_map_2, 'completness', comp_2, 'selective', [1 5 3 7])
%
% 4. neighbor
% vol_viz(gid_map_1, 'neighbor')
% vol_viz(gid_map_2, 'neighbor','selective', [1 5 3 7])
%
% 5. exclude & selective
% vol_viz(gid_map_2, 'orientation', grain_rodV_2, 'cubic', 'selective', grain_surface_2)
% vol_viz(gid_map_2, 'neighbor', 'exclude', grain_surface_2)
%==========================================================================

gid_list = unique(gid_map(:));
gid_list = gid_list(2:end);
color_gid_list = zeros(length(gid_list),3);

%% define selective grains
if any(strcmp(varargin,'selective')) 
    idx = find(strcmp(varargin,'selective'))+1;
    selective = varargin{idx};
elseif any(strcmp(varargin,'exclude')) 
    idx = find(strcmp(varargin,'exclude'))+1;
    exclude = varargin{idx};
end

if exist('selective','var')
    gid_list(~ismember(gid_list,selective))=[];
    no_grain = selective(~ismember(selective,gid_list));
    no_grain = reshape(no_grain,[1,length(no_grain)]);
    fprintf('lisf of non-exsiting grain ID from input grain ID.\n');
    disp(no_grain);
elseif exist('exclude','var')
    no_grain = exclude(~ismember(exclude,gid_list));
    no_grain = reshape(no_grain,[1,length(no_grain)]);
    gid_list(ismember(gid_list,exclude))=[];
    fprintf('lisf of non-exsiting grain ID from input grain ID.\n');
    disp(no_grain);
end

%% define coloring

if ~any(strcmp(varargin,'orientation') | ...
        strcmp(varargin,'completness') | ...
        strcmp(varargin,'size') | ...
        strcmp(varargin,'neighbor'))
    %random coloring
    color_gid_list = rand(length(gid_list),3);
        
    
elseif any(strcmp(varargin,'orientation')) 
    %coloring based on crystallographic direction
    idx = find(strcmp(varargin,'orientation'))+1;
    grain_rodV = varargin{idx};
    crystal = varargin{idx+1};
    cs = crystalSymmetry(crystal);
    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;
        if exist('selective','var')
            grain_rodV(~ismember(gid_list,selective))=[];
        elseif exist('exclude','var')
            grain_rodV(ismember(gid_list,exclude))=[];
        end
    for i = 1:length(gid_list)
        q1 = rodrigues2quat(vector3d(grain_rodV(i,:)));
        o1 = orientation(q1, cs);
        if isnan(o1.phi1)
            color_gid_list(i,:) = [0 0 0];
        else
            color_gid_list(i,:) = oM.orientation2color(o1);
        end
    end
    
    
elseif any(strcmp(varargin,'completness')) 
    %coloring based on average completeness
    idx = find(strcmp(varargin,'completness'))+1;
    comp = varargin{idx};
        if exist('selective','var')
            comp(~ismember(gid_list,selective),:)=[];
        elseif exist('exclude','var')
            comp(ismember(gid_list,exclude))=[];
        end
    comp_mean = zeros(length(gid_list),1);
    for i = 1:length(gid_list)
        comp_mean(i) = mean(comp(gid_map==gid_list(i)));
    end
    comp_color = jet(101);
    for i = 1:length(gid_list)
        color_gid_list(i,:) = ...
            comp_color(ceil( ( comp_mean(i)-min(comp_mean)+0.000001 ) / ( max(comp_mean) -min(comp_mean) ) *100 ),:);
    end
    
    
elseif any(strcmp(varargin,'size')) 
    %coloring based on grain size
    idx = find(strcmp(varargin,'size'))+1;
    numElement = varargin{idx};
        if exist('selective','var')
            numElement(~ismember(gid_list,selective),:)=[];
        elseif exist('exclude','var')
            numElement(ismember(gid_list,exclude))=[];
        end
    size_color = jet(max(numElement(:,2)));
    for i = 1:length(gid_list)
        color_gid_list(i,:) = size_color(numElement(i,2),:);
    end

elseif any(strcmp(varargin,'neighbor')) 
    %coloring based on the number of neighbor grain
    adj = imRAG(gid_map);
    adj((adj(:,1) > adj(:,2)),:) = [];
    [~,F] = mode(adj(:));
    neighbor_color = jet(F);
    for i = 1:length(gid_list)
        numNeighbor = sum(sum(adj == gid_list(i)));
        color_gid_list(i,:) = neighbor_color(numNeighbor,:);
    end
end

%% visualization

figure
for i = 1:length(gid_list)
    fprintf('Visualizing %d out of %d grains\n', i, numel(gid_list));
    FV = isosurface(gid_map == gid_list(i));
    
    if size(FV.vertices,1) > 0
        patch(smoothpatch(FV,1,3), ...
                'facecolor', color_gid_list(i,:), ...
                'edgecolor', 'none', ...
                'facealpha', 1);
    end
        hold on;
end


view([30,20]); axis equal; camlight left; camlight right; material metal;

%different colorbar for different coloring
    if any(strcmp(varargin,'completness'))
        colormap(comp_color)
        caxis([min(comp_mean) max(comp_mean)])
        colorbar
    elseif any(strcmp(varargin,'size'))
        colormap(size_color)
        colorbar
        caxis([0 max(numElement(:,2))])
    elseif any(strcmp(varargin,'neighbor'))    
        colormap(neighbor_color)
        colorbar
        caxis([0 F]) 
    end
    
end
