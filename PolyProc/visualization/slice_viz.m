function [] = slice_viz(gid_map, varargin)
% slice_visualization reconstruct specified 2d slice
%==========================================================================
% FILENAME:          slice_visualization.m
% DATE:              1 May, 2019        
% PURPOSE:           reconstruct 2d slice
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%    z_pos     : slice of interest (along sample's axis of rotation)
%
%OPTIONAL : specify grain coloring criterion 
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
%OUT :
%    2D reconstruction of a slice
%==========================================================================
%EXAMPLE :
% 0. no specification
%    slice_viz(gid_map_1,'slice',77)
%    slice_viz(gid_map_1,'slice',[70 71 72 73 74 75])
%
% 1. orientation
%    slice_viz(gid_map_1,'slice',77,'orientation', grain_rodV_1, 'cubic')
%    slice_viz(gid_map_1,'slice',[70 71 72 73 74 75],'orientation', grain_rodV_1, 'cubic')
%
% 2. size
%    slice_viz(gid_map_1,'slice',77,'size', numElement_1)
%    slice_viz(gid_map_1,'slice',[70 71 72 73 74 75],'size', numElement_1)
%
% 3. completeness
%    slice_viz(gid_map_1,'slice',77,'completness', comp_1)
%
% 4. neighbor
%    slice_viz(gid_map_1,'slice',77,'neighbor')
%
%==========================================================================

%% set up different z_pos
if any(strcmp(varargin,'slice')) 
    idx = find(strcmp(varargin,'slice'))+1;
    z_pos = varargin{idx};
else 
    z_pos = varargin{1};
end

%% set up background

gid_list_tot = unique(gid_map(:));
gid_list_tot = gid_list_tot(2:end);

gid_slice = zeros(size(gid_map,1),size(gid_map,2),length(z_pos));
for j = 1:length(z_pos)
    gid_slice(:,:,j) = gid_map(:,:,z_pos(j));
end

    gid_slice_or = zeros(size(gid_slice,1), size(gid_slice,2), 3);
    gid_slice_tot = zeros(size(gid_slice_or));

gid_list = unique(gid_slice(:));
gid_list = gid_list(2:end);
color_gid_list = zeros(size(gid_list,1), 3);


%% define coloring

if ~any(strcmp(varargin,'orientation') | ...
        strcmp(varargin,'completness') | ...
        strcmp(varargin,'size') | ...
        strcmp(varargin,'neighbor'))
    %random coloring
    color_gid_list = rand(size(gid_list,1),3);
    
elseif any(strcmp(varargin,'orientation')) 
    %coloring based on crystallographic direction
    idx = find(strcmp(varargin,'orientation'))+1;
    grain_rodV = varargin{idx};
    crystal = varargin{idx+1};
    cs = crystalSymmetry(crystal);
    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;
    for i = 1:length(gid_list)
        q1 = rodrigues2quat(vector3d(grain_rodV((gid_list_tot==gid_list(i)),:)));
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
    comp_mean = zeros(length(gid_list),1);
    for i = 1:length(gid_list)
        comp_mean(i) = mean(comp(gid_map==gid_list(i)));
    end
    comp_color = jet(101);
    for i = 1:length(gid_list)
        color_gid_list(i,:) = ...
            comp_color(ceil( ( comp_mean(i)-min(comp_mean)+0.001 ) / ( max(comp_mean) -min(comp_mean) ) *100 ),:);
    end
    
    
elseif any(strcmp(varargin,'size')) 
    %coloring based on grain size
    idx = find(strcmp(varargin,'size'))+1;
    numElement = varargin{idx};
    size_color = jet(max(numElement(:,2)));
    for i = 1:length(gid_list)
        color_gid_list(i,:) = size_color(numElement(gid_list_tot==gid_list(i),2),:);
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

for k = 1:length(z_pos)
            gid_slice_layer = gid_slice(:,:,k);
            %change background color to white
            testIm = gid_slice_layer == 0;
            gid_slice_or_copy = gid_slice_or;
            gid_slice_or_copy(:,:,1) =  testIm;
            gid_slice_or_copy(:,:,2) =  testIm;
            gid_slice_or_copy(:,:,3) =  testIm;
            gid_slice_tot_copy = gid_slice_tot + gid_slice_or_copy;


            fprintf('Displaying slice #%d ...\n', z_pos(k));
            gid_list_layer = unique(gid_slice_layer(:));
            gid_list_layer = gid_list_layer(2:end);
            
            for i = 1:length(gid_list_layer)

                testIm = (gid_slice_layer == gid_list_layer(i));
                
                gid_slice_or_copy(:,:,1) = color_gid_list(gid_list==gid_list_layer(i),1) .* testIm;
                gid_slice_or_copy(:,:,2) = color_gid_list(gid_list==gid_list_layer(i),2) .* testIm;
                gid_slice_or_copy(:,:,3) = color_gid_list(gid_list==gid_list_layer(i),3) .* testIm;

                gid_slice_tot_copy = gid_slice_tot_copy + gid_slice_or_copy;
            end

            figure, imagesc(gid_slice_tot_copy);
            axis equal;
            axis off;
            % colorbar setup
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


end
