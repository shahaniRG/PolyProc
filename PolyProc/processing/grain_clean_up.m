function [gid_map, rodV, comp, numElement, grain_rodV, grain_coord, grain_surface, adj, cs] = ...
                    grain_clean_up(file_name,varargin)
% grain_clean_up process DCT data set 
% to cluster subgrains into grains, 
% get rid of noise voxels, 
% and exclude untrustworthy grain.     
%==========================================================================
% FILENAME:          grain_clean_up.m
% DATE:              1 May, 2019     
% PURPOSE:           process DCT data set
%==========================================================================
%IN :
%    file_name  : (string) filename of hdf5 file (.h5)
%
%OPTIONAL0 : (string) specify the crystal structure for crystallographic analysis
%     crystal
%     |---'cubic' (default)
%     |---further_update  - to be updated
%OPTIONAL1 : in case of the inputting one of multiple time step data, 
%            and registration is needed
%    mask        : (logical, array) 3-dimensional data set, determining volume of interest
%                   default : no mask, all true
%    tform       : 
%     |---(1*1 affine3d) contains transformation matrix (voxel tranformation)
%          default : no registration on scalar data(no translation, no rotation)
%     |---(1*1 vector3d) contains rotation rodrigues vectors of each grain direction
%          default : no rotation on vector data
%
%    misoriThresh       : (double) angular threshold
%                         default : 1
%    grain_size_mini    : (double) grain size threshold
%                         default : 1
%    completeness_mini  : (double) completeness threshold (value btw 0 & 1)
%                         default : 0.4
%OUT : 
%    gid_map    : (array of 3D dataset) cleaned up gid_map
%    rodV       : (array of 4D dataset) cleaned up rodrigues vectors
%    comp       : (array of 3D dataset) cleaned up completeness
%    numElement : n*2 array with (1st coloumn-gid) & (2nd column-numbner of volexs)
%    grain_rodV : n*3 array with rodrigues vectors of grains
%    grain_coord: n*3 array with center of mass coordination
%    grain_surface: n*1 array specifying surface touching grains
%    adj        : n*2 array shows grain adjacency in number ascending order
%    cs         : (crystalSymmetry) cell containing crystal structrue & symmetry information
%==========================================================================
%EXAMPLE 1 : dealing with multiple time steps, so tform & mask is ready
%           [gid_map, rodV, comp, numElement, grain_rodV, grain_coord, adj, cs]...
%              grain_clean_up('t2_1.h5',...
%              'crystal','cubic',...
%              'tform', tform_matrix_t2, rot_rod_t2, ...
%              'mask',gid_map_mask, ...
%              'threshold_ang',1,'threshold_vol',33,'threshold_comp',0.1);
%
%EXAMPLE 2 : dealing with single time step
%           [gid_map, rodV, comp, numElement, grain_rodV, grain_coord, adj, cs]...
%              grain_clean_up('t2_1.h5','threshold_ang',1,'threshold_vol',33,'threshold_comp',0.1);
%
%EXAMPLE 3 : with minimal input variable
%           [gid_map, rodV, comp, numElement, grain_rodV, grain_coord, adj, cs]...
%              grain_clean_up('t2_1.h5');
%
%==========================================================================


    % Read in Grain Id Map, rodrigues vectors, voxel size, and completeness
    gid_map = h5read(file_name,'/LabDCT/Data/GrainId');
    rodV = h5read(file_name,'/LabDCT/Data/Rodrigues');
    comp = double(h5read(file_name,'/LabDCT/Data/Completeness'));
    
   %% Default values
   defaultCrystal = 'cubic';
   expectedCrystal = {'cubic','further_update'};
   defaultmask = zeros(size(gid_map))+1;
   defaultTform = affine3d([1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ;  0 0 0 1 ]);
   defaultRot = Rodrigues(rotation('axis',xvector,'angle',0*degree));
   defaultMisoriThresh = 1;
   defaultGrain_size_mini = 1;
   defaultCompleteness_mini = 0.1;
   
   %% Define Crystal Structure
    if any(strcmp(varargin,'crystal'))
        idx = find(strcmp(varargin,'crystal'))+1;
        crystal= expectedCrystal{strcmp(expectedCrystal, varargin{idx})};
    else
        crystal = defaultCrystal;
    end
    
    % define crystal structure & symmetry
    cs = crystalSymmetry(crystal);
    
    %% Define Transformation Matrix
    if any(strcmp(varargin,'tform'))
        idx = find(strcmp(varargin,'tform'))+1;
        tform = varargin{idx};
        rot_rod = varargin{idx+1};
    else
        tform = defaultTform;
        rot_rod = defaultRot;
    end
    
    %% Define Mask Matrix
    if any(strcmp(varargin,'mask'))
        idx = find(strcmp(varargin,'mask'))+1;
        mask = varargin{idx};
    else
        mask = defaultmask;
    end
    
    %% Define Threshold Values
    if any(strcmp(varargin,'threshold_ang'))
        idx = find(strcmp(varargin,'threshold_ang'))+1;
        misoriThresh = varargin{idx};
    else
        misoriThresh = defaultMisoriThresh;
    end
    
    if any(strcmp(varargin,'threshold_vol'))
        idx = find(strcmp(varargin,'threshold_vol'))+1;
        grain_size_mini = varargin{idx};
    else
        grain_size_mini = defaultGrain_size_mini;
    end
    
    if any(strcmp(varargin,'threshold_comp'))
        idx = find(strcmp(varargin,'threshold_comp'))+1;
        completeness_mini = varargin{idx};
    else
        completeness_mini = defaultCompleteness_mini;
    end

    fprintf('Thresholds are identified.\n')
    %% Applying Transformation (translation + roatation) + Rodrigues update + Mask Matrix = Setting Region of Interest

    % Transformation (translation + rotation)
    
    % (1) apply scope mask to grain ID
    gid_map = imwarp(gid_map,tform,'OutputView',imref3d(size(mask)),'interp','nearest');
    gid_map = double(gid_map).*mask;
    
    % (2) apply scope mask to rodrigues vectors
    r1 = squeeze(rodV(1,:,:,:));
    r2 = squeeze(rodV(2,:,:,:));
    r3 = squeeze(rodV(3,:,:,:));
    r1 = imwarp(r1,tform,'OutputView',imref3d(size(mask)),'interp','nearest');
    r1 = double(r1).*mask;
    r2 = imwarp(r2,tform,'OutputView',imref3d(size(mask)),'interp','nearest');
    r2 = double(r2).*mask;
    r3 = imwarp(r3,tform,'OutputView',imref3d(size(mask)),'interp','nearest');
    r3 = double(r3).*mask;
    rodV = permute(cat(4,r1,r2,r3),[4,1,2,3]);
    
    clear r1 r2 r3
    
    % (3) apply scope mask to completeness
    comp = imwarp(comp,tform,'OutputView',imref3d(size(mask)),'interp','nearest');
    comp = double(comp).*mask;    

    
    % Rodrigues update
    
    % (4) apply rotation to each rodrigues elements
    rodV(1,:,:,:) = rodV(1,:,:,:) + rot_rod.x;
    rodV(2,:,:,:) = rodV(2,:,:,:) + rot_rod.y;
    rodV(3,:,:,:) = rodV(3,:,:,:) + rot_rod.z;
    
    fprintf('Data have been transformed.\n')
    
    %% detect (surface + top + bottom) grains on exterior
    gid_map_copy = gid_map;
    gid_map_filled = imfill(logical(gid_map_copy),'holes');
    gid_map_copy(gid_map_filled==0) = -100;
    
    adj_interior = imRAG(gid_map_copy);
    grain_surface = adj_interior(adj_interior(:,1)==-100,2);   
    
    fprintf('Surface grains are identified.\n')
    %% Clean up by Angular Threshold
    fprintf('Clustering data by angular threshold . . .\n')
    
    % Compute adjacency matrix
    adj = imRAG(gid_map);
    adj((adj(:,1) > adj(:,2)),:) = [];
    unique_gid = unique(adj);
   
    % Get grain-averaged Rodrigues vectors
    grain_rodV = zeros(max(unique_gid),3);
    for i = transpose(unique_gid)
        [x,y,z] = ind2sub(size(gid_map),find(gid_map == i));
        grain_rodV(i,1) = nanmean(rodV(sub2ind(size(rodV), ...
                    1*ones(size(x)),x,y,z)));
        grain_rodV(i,2) = nanmean(rodV(sub2ind(size(rodV), ...
                    2*ones(size(x)),x,y,z)));
        grain_rodV(i,3) = nanmean(rodV(sub2ind(size(rodV), ...
                    3*ones(size(x)),x,y,z)));
    end
    
    % Initialize some variables
    gid_map_new = gid_map;
    noGrains = numel(unique(gid_map_new(:))); % this contains '0'(space in the unique grain list)
    iter = 1;
    goodPairs = [0,0]; % Array that holds grain i,j pairs 
                       % that have misorientation > misoriThresh

    % Now the fun stuff: clustering operation
    while iter < size(adj,1)

        % Pick two touching grains in adj set
        i = adj(iter,1); % adj(quickInd,1);
        j = adj(iter,2); % adj(quickInd,2); 

        if ~ismember([i,j], goodPairs, 'rows')

            fprintf(' No. Grains = %d. Checking Grains %d and %d ... \n', ...
                    noGrains-1, i, j);

            % Compute their misorientation
            q1 = rodrigues2quat(vector3d(grain_rodV(i,:)));
            o1 = orientation(q1, cs);
            q2 = rodrigues2quat(vector3d(grain_rodV(j,:)));
            o2 = orientation(q2, cs);
            misori = angle(o1,o2)/degree;

            % If their misorientation is less than pre-set value, then ...
            if misori < misoriThresh

                % (i) Combine the smaller grain into the larger 
                if sum(gid_map_new(:) == i) > sum(gid_map_new(:) == j)
                   gid_map_new(gid_map_new == j) = i;
                else
                   gid_map_new(gid_map_new == i) = j;
                end
                iter = 1;

                % (ii) Recompute the number of grains after clustering
                noGrains = noGrains - 1;

                % (iii) Recompute the grains that touch (adjacent grains)
                adj = imRAG(gid_map_new);
                adj((adj(:,1) > adj(:,2)),:) = [];

            else
                % Keep marching through the list of neighbors in adj
                iter = iter+1;
                goodPairs = [goodPairs; i,j];

            end

        else
            % Keep marching through the list of neighbors in adj
            iter = iter+1;
        end 
    end
    fprintf('Angular thresholding is complete. Number of grains = %d\n', noGrains-1);
    gid_map = gid_map_new;
    % calculate volume after angular thershold
    % number_voxel1 = length(find(gid_map~=0));
    clear x y z gid_map_new noGrains iter goodPairs i j q1 o1 q2 o2 misori misoriThresh
      
    %% Clean up by Volume Threshold
    
    gid_mask = gid_map;
    gid_mask(gid_mask == 0) = [];
    numElement = accumarray(gid_mask(:),1);
    numElement_in = [1:length(numElement)]';
    numElement = [numElement_in, numElement];
    numElement_in(numElement(:,2) >= grain_size_mini) = [];
    gid_map(ismember(gid_map,numElement_in)) = 0;
    
    % update adj and numElement
    numElement(numElement(:,2) < grain_size_mini,:) = [];
    adj = imRAG(gid_map);
    adj((adj(:,1) > adj(:,2)),:) = [];
    unique_gid = unique(adj);
    numElement(~ismember(numElement(:,1),unique_gid),:) = [];
    
    clear gid_mask numElement_in grain_size_mini
    fprintf('Volume thresholding is complete.\n')
    
    %% Clean up by Completeness Threshold
    
    comp_mean = zeros(length(numElement),1);
    % calculate completeness of grain (voxel average completeness)
    for i = 1:length(numElement(:,1))
        comp_mean(i) = mean(comp(gid_map==numElement(i,1)));
    end
    % remove grains with low completeness
    untrust_grain = numElement((0<comp_mean & comp_mean<completeness_mini),1);
    gid_map(ismember(gid_map,untrust_grain))=0;
    numElement(ismember(numElement(:,1),untrust_grain),:)=[];
    
    % update adj
    adj = imRAG(gid_map);
    adj((adj(:,1) > adj(:,2)),:) = [];
    
    clear comp_mean i completeness_mini untrust_grain
    fprintf('Completeness thresholding is complete.\n')

    %% Update variables
    grain_rodV = zeros(length(numElement(:,1)),3);
    grain_coord = grain_rodV;
    % update grain-averaged Rodrigues vectors    
    for i = 1:length(numElement(:,1))
        [x,y,z] = ind2sub(size(gid_map),find(gid_map==numElement(i,1)));
        grain_coord(i,:) = round(mean([x,y,z]));
        grain_rodV(i,1) = nanmean(rodV(sub2ind(size(rodV), ...
                1*ones(size(x)),x,y,z)));
        grain_rodV(i,2) = nanmean(rodV(sub2ind(size(rodV), ...
                2*ones(size(x)),x,y,z)));
        grain_rodV(i,3) = nanmean(rodV(sub2ind(size(rodV), ...
                3*ones(size(x)),x,y,z)));
    end
    
end