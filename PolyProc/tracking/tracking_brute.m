function [grain_matching, match_eff] = tracking_brute(...
                            gid_map_1, numElement_1, grain_rodV_1,...
                            gid_map_2, numElement_2, grain_rodV_2,...
                            varargin)
% tracking_brute tracks grains between two different time steps
%==========================================================================
% FILENAME:          tracking_brute.m
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
%    [grain_matching_brute, matching_efficiency_brute] = ...
%                     tracking_brute(gid_map_1, ...
%                     numElement_1, grain_rodV_1, ...
%                     gid_map_2, numElement_2, grain_rodV_2,...
%                     'cubic', 2, 4);   
%==========================================================================                     
                        
tic
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
%% 
    % setting up variables
    grain_matching = [numElement_1(:,1),zeros(size(numElement_1,1),1)];
    candidate_storage = zeros(length(numElement_1),200);
    misori_tot = candidate_storage+100;
    dist_centroid_tot = misori_tot;
    cs = crystalSymmetry(crystal);
%     gid_map_2_original = gid_map_2;
    
    for j = 1:length(numElement_1(:,1))
        fprintf('Matching grain %d out of %d total grains ... \n', ...
                                j, length(numElement_1(:,1)));
        [x,y,z] = ind2sub(size(gid_map_1),find(gid_map_1 == numElement_1(j,1)));
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

        % define centroid and orientation of grain at 'n' step
        centroid_1 = round(mean([x,y,z]));
        q1 = rodrigues2quat(vector3d(grain_rodV_1(j,:)));
        o1 = orientation(q1, cs);
        
        %apply on 'n+1' step
        gid_scope = area_mask_1.*gid_map_2;
        candidate = unique(gid_scope);
        candidate(candidate == 0) = [];
        candidate_storage(j,1:length(candidate))=candidate;

        % define centroid and orienetation of grina at 'n+1' step
        for k = 1:length(candidate)
            ind = numElement_2(:,1)==candidate(k);
            
            [x2,y2,z2] = ind2sub(size(gid_map_2),find(gid_map_2 == candidate(k)));
            centroid_2 = round(mean([x2,y2,z2]));
            dist_centroid = sqrt(sum(((centroid_1-centroid_2).^2),2));
            dist_centroid_tot(j,k) = dist_centroid;
            
            q2 = rodrigues2quat(vector3d(grain_rodV_2(ind,:)));
            o2 = orientation(q2, cs);
            misori = angle(o2,o1)/degree;
            misori_tot(j,k) = misori;
        end

        ind_ang = find(misori_tot(j,:)< misori_threshold);
        ind_dist = find(dist_centroid_tot(j,:) < dist_threshold);
        
        if isempty(intersect(ind_ang,ind_dist))
            grain_matching(j,2) = NaN;
        else
            ind_candidate_final=intersect(ind_ang,ind_dist);
            cost_ang = misori_tot(j,ind_candidate_final);
            cost_dist = dist_centroid_tot(j,ind_candidate_final);
            cost = (1-weight_factor).*cost_dist + (weight_factor).*cost_ang;
            ind_cost = find(cost==min(cost));
            grain_matching(j,2) = candidate_storage(j,ind_candidate_final(ind_cost));
            gid_map_2(gid_map_2==candidate_storage(j,ind_candidate_final(ind_cost)))=0;
        end

    end
    
    grain_matching_2 = grain_matching(~isnan(grain_matching(:,2)),:);
    match_eff = size(grain_matching_2,1)/length(numElement_2);
    toc 
    fprintf(' Assignments made with efficiency %d.\n', match_eff); 

end

