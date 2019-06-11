function [grain_matching, match_eff] = ...
                            tracking_hun(...
                            gid_map_1, numElement_1, grain_rodV_1,...
                            gid_map_2, numElement_2, grain_rodV_2,...
                            varargin)
% grain_tracking_hungarian tracks grains between two different time steps
%==========================================================================
% FILENAME:          grain_tracking_hungarian.m
% DATE:              1 May, 2019        
% PURPOSE:           grain tracking
%==========================================================================
%IN :
%    gid_map_1   : (array) previous time step 3D data set of gid_map
%
%    gid_map_2   : (array) after time step 3D data set of gid_map
%
%    numElement_1&2 : n*2 array with (1st coloumn-gid) & (2nd column-numbner of voxels)
%
%    grain_rodV_1&2 : n*3 array with rodrigues vectors of grains
%
%OPTIONAL 1: 
%    crystal      : (string) crystal stucture of systems
%                   (default - 'cubic')
%
%    padding   : (double) number of voxels to extend for neighboring grain investigation
%                   (distance theshold is automatically defined based on volume of each grain & area_scope value)
%                   (default - 2)
%
%    misori_threshold : (double) misorientation threshold to be considered as potential matching grain
%                       (default - 5)
%
%OPTIONAL 2: 
%    weight_factor : (double) value between 0 and 1.
%                    weight factor of misorientation over distance.
%                    cost matrix is calculated based on the weight factor.
%                    (1-weight_factor automatically becomes weight_factor of distance)
%                    (default : 1)
%    visualization flag : (string) 'viz' with number of grains to visualize.
%                         (double) the following integer represent the
%                                  number of grain to visualize
%                                    (randomly display matched grains)
%                         (default : no visualization)
%
%OUT :
%    grain_matching : n*2 array with grain paring (1st coloumn-previous time step gid)
%                                                 (2nd coloumn-next time step gid)
%
%    match_eff      : (double) matching efficiency of grain tracking
%                      = (number of paired grains in next time step object) / 
%                         (entire number of grains in next time step object)
%
%==========================================================================
%EXAMPLE :
%    [grain_matching_hun_1_opt, matching_efficiency_hun_1_opt] = ...
%                     tracking_hun...
%                     (gid_map_1, numElement_1, grain_coord_1, grain_rodV_1,...
%                     gid_map_2, numElement_2, grain_coord_2, grain_rodV_2,...
%                     'cubic', 'extension', 2, 'miorientation', 10);
%==========================================================================
tic

fprintf('Defining grain center-of-mass.\n');
% define center of mass for every grain (for later distance calculation)
grain_coord_1 = zeros(length(numElement_1(:,1)),3);
grain_coord_2 = zeros(length(numElement_2(:,1)),3);
    for i = 1:length(numElement_1(:,1))
        [x,y,z] = ind2sub(size(gid_map_1),find(gid_map_1==numElement_1(i,1)));
        grain_coord_1(i,:) = round(mean([x,y,z]));
    end
    
    for j = 1:length(numElement_2(:,1))
        [x,y,z] = ind2sub(size(gid_map_2),find(gid_map_2==numElement_2(j,1)));
        grain_coord_2(j,:) = round(mean([x,y,z]));
    end

    %% optional input1 - crystal, area_scope, misori_threshold
    
    if any(strcmp(varargin,'crystal'))
        idx = find(strcmp(varargin,'crystal'))+1;
        crystal= varargin{idx};
    else
        crystal = 'cubic';
    end
    
    if any(strcmp(varargin,'extension'))
        idx = find(strcmp(varargin,'extension'))+1;
        area_scope= varargin{idx};
    else
        area_scope = 2;
    end
    
    if any(strcmp(varargin,'misorientation'))
        idx = find(strcmp(varargin,'misorientation'))+1;
        misori_threshold= varargin{idx};
    else
        misori_threshold = 5;
    end
    
    
    %% optional input2 - weight factor between distance and misorientation
    if any(strcmp(varargin,'weight_factor'))
        idx = find(strcmp(varargin,'weight_factor'))+1;
        weight_factor= varargin{idx};
    else
        weight_factor = 1;
    end
    
    if any(strcmp(varargin,'viz'))
        idx = find(strcmp(varargin,'viz'))+1;
        vis_grain= varargin{idx};
        flag = 1;
    else
        flag = 0;
    end
    
    
%% 
    cs = crystalSymmetry(crystal);
    
    fprintf('Setting up cost matrix ...\n');
    % setting rectangular cost matrices for hungarian algorithm
    cost_mat = Inf(length(numElement_1),length(numElement_2));
    cost_mat_dist = cost_mat;
    cost_mat_ori = cost_mat;
    
    
    for i = 1:length(numElement_1)
        
                %setting up the scope
                [x,y,z] = ind2sub(size(gid_map_1),find(gid_map_1 == numElement_1(i,1)));
                x_min = (min(x)-area_scope);
                x_max = (max(x)+area_scope);
                y_min = (min(y)-area_scope);
                y_max = (max(y)+area_scope);
                z_min = (min(z)-area_scope);
                z_max = (max(z)+area_scope);
                    if x_min < 1
                        x_min = 1;
                    end
                    if x_max > size(gid_map_1,1)
                        x_max = size(gid_map_1,1);
                    end
                    if y_min < 1
                        y_min = 1;
                    end
                    if y_max > size(gid_map_1,2)
                        y_max = size(gid_map_1,2);
                    end
                    if z_min < 1
                        z_min = 1;
                    end
                    if z_max > size(gid_map_1,3)
                        z_max = size(gid_map_1,3);
                    end
                area_mask_1 = zeros(size(gid_map_1));    
                area_mask_1( x_min:x_max, y_min:y_max, z_min:z_max) = 1;
                dist_threshold = 0.5*sqrt(sum((([x_max, y_max, z_max]-[x_min,y_min,z_min]).^2),2));
                gid_scope = area_mask_1.*gid_map_2;
                candidate = unique(gid_scope);
                candidate(candidate == 0) = [];
        
        
        q1 = rodrigues2quat(vector3d(grain_rodV_1(i,:)));
        o1 = orientation(q1, cs);
        
        for j = candidate'
            
            if ~isempty(candidate)
            
            % calculate distance btw grain centroid
            dist = sqrt( sum( ( grain_coord_1(i,:) - grain_coord_2(numElement_2(:,1)==j,:) ).^2 , 2 ) );
            if dist < dist_threshold

                q2 = rodrigues2quat(vector3d(grain_rodV_2(numElement_2(:,1)==j,:)));
                o2 = orientation(q2, cs);
                ori = angle(o2,o1)/degree;

                if dist < dist_threshold && ori < misori_threshold
                cost_mat_dist(i,numElement_2(:,1)==j) = dist/dist_threshold;
                cost_mat_ori(i,numElement_2(:,1)==j) = ori/misori_threshold;
                end

            end
            end
        
        end
    end
    
    cost_mat = (1-weight_factor)*cost_mat_dist + (weight_factor)*cost_mat_ori;
    
    fprintf('Performing combinatorial optimization.\n');
    [assignment, ~] = assignmentsuboptimal1(cost_mat);
    
    grain_matching_2 = [numElement_1(:,1) assignment];

    for i = 1:length(grain_matching_2)
        if assignment(i) ~= 0
            grain_matching_2(i,2) = numElement_2(assignment(i),1);
        else
            grain_matching_2(i,2) = NaN;
        end
    end
    
    grain_matching = grain_matching_2;
toc
    
    grain_matching_2 = grain_matching(~isnan(grain_matching(:,2)),:);
    match_eff = size(grain_matching_2,1)/length(numElement_2);
    fprintf(' Assignments made with efficiency %d.\n', match_eff);
%%    grain match visualization
if flag == 1
    fprintf('Visualizing matching grains ...\n');
    figure
    subplot(1,2,1)
    S_A = isosurface(gid_map_1,0);
    patch(S_A, 'Facecolor', 'b', 'EdgeColor', 'none','FaceAlpha',0.1);
    subplot(1,2,2)
    S_C = isosurface(gid_map_2,0);
    patch(S_C, 'Facecolor', 'b', 'EdgeColor', 'none','FaceAlpha',0.1);

    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;

    color_gid_list_1 = zeros(length(grain_matching_2), 3);
    color_gid_list_2 = color_gid_list_1;


    for i =  randi([1 size(grain_matching_2,1)],1,vis_grain)

            idx_1 = numElement_1(:,1)==grain_matching_2(i,1);
            q1 = rodrigues2quat(vector3d(grain_rodV_1(idx_1,:)));
            o1 = orientation(q1, cs);
            if isnan(o1.phi1)
                color_gid_list_1(i,:) = [0 0 0];
            else
                color_gid_list_1(i,:) = oM.orientation2color(o1);
            end
            subplot(1,2,1)
            hold on
            S_B = isosurface(gid_map_1==grain_matching_2(i,1),0);
            if size(S_B.vertices,1) > 0
                patch(smoothpatch(S_B,1,3), 'Facecolor', color_gid_list_1(i,:),...
                        'EdgeColor', 'none','FaceAlpha',0.8);
            end

            idx_2 = numElement_2(:,1)==grain_matching_2(i,2);
            q2 = rodrigues2quat(vector3d(grain_rodV_2(idx_2,:)));
            o2 = orientation(q2, cs);
            if isnan(o2.phi1)
                color_gid_list_2(i,:) = [0 0 0];
            else
                color_gid_list_2(i,:) = oM.orientation2color(o2);
            end
            subplot(1,2,2)
            hold on
            S_D = isosurface(gid_map_2==grain_matching_2(i,2),0);
            if size(S_B.vertices,1) > 0
            patch(smoothpatch(S_D,1,3), 'Facecolor', color_gid_list_2(i,:), ...
                    'EdgeColor', 'none','FaceAlpha',0.8);
            end
    end
        subplot(1,2,1)
        title('t_{1}');
        view([30,20]); axis equal; camlight left; camlight right; material metal;
        axis off
        subplot(1,2,2)
        title('t_{2}');
        view([30,20]); axis equal; camlight left; camlight right; material metal;
        axis off

end

end

