function [grain_matching, match_eff] = ...
                            tracking_hun(...
                            numElement_1, grain_rodV_1,gcenterArray1,...
                            numElement_2, grain_rodV_2,gcenterArray2,...
                            CandidateCArray, dist_thresh_array,...
                            varargin)
%% Optional input1 - crystal, area_scope, misori_threshold
    
    if any(strcmp(varargin,'crystal'))
        idx = find(strcmp(varargin,'crystal'))+1;
        crystal= varargin{idx};
    else
        crystal = 'cubic';
    end
    
    if any(strcmp(varargin,'misorientation'))
        idx = find(strcmp(varargin,'misorientation'))+1;
        misori_threshold= varargin{idx};
    else
        misori_threshold = 5;
    end

%% Tracking with Cost Matrix

    cs = crystalSymmetry(crystal);
    fprintf('Setting up cost matrix ...\n');

    %Setting rectangular cost matrices for hungarian algorithm
    cost_mat      = Inf(length(numElement_1),length(numElement_2));
    cost_mat_dist = cost_mat;
    cost_mat_ori  = cost_mat;
    
    for i = 1:length(numElement_1)
        fprintf('Grain %d of %d...\n', i, length(numElement_1));
        candidate      = CandidateCArray{i};
        dist_threshold = 2*dist_thresh_array(i);
        gIDcheck       = numElement_1(i,1);

        %Orientation for initial grain
        q1 = rodrigues2quat(vector3d(grain_rodV_1(i,:)));
        o1 = orientation(q1, cs);
        
        if ~isempty(candidate)
            for j = candidate'  
                
                %Calculate distance btw grain centroid
                dist = sqrt((gcenterArray1(gIDcheck,1) - gcenterArray2(j,1))^2 +...
                            (gcenterArray1(gIDcheck,2) - gcenterArray2(j,2))^2 +...
                            (gcenterArray1(gIDcheck,3) - gcenterArray2(j,3))^2);
                
                %Calculate volume overlap
                %overlap_frac = overlap_matrix(j,gIDcheck);

                %Calculate Misori  
                if ~isempty(grain_rodV_2(numElement_2(:,1)==j,:))
                    q2  = rodrigues2quat(vector3d(grain_rodV_2(numElement_2(:,1)==j,:)));
                    o2  = orientation(q2, cs);
                    ori = angle(o1,o2)/degree;
                    
                    %Populate Cost Matrix
                    if ori<misori_threshold
                        cost_mat_dist(i,numElement_2(:,1)==j)  = dist/dist_threshold;
                        cost_mat_ori(i,numElement_2(:,1)==j)   = ori/misori_threshold;
                    end
                end
            end        
        end
    end

    weight_factor_dist = 0.5;
    cost_mat = (weight_factor_dist)*cost_mat_dist + (1-weight_factor_dist)*cost_mat_ori;
    
    fprintf('Cost Matrix Populated - Performing combinatorial optimization.\n');
    [assignment, ~] = assignmentsuboptimal1(cost_mat);
    
    grain_matching_temp = [numElement_1(:,1) assignment];

    for i = 1:length(grain_matching_temp)
        if assignment(i) ~= 0
            grain_matching_temp(i,2) = numElement_2(assignment(i),1);
        else
            grain_matching_temp(i,2) = NaN;
        end
    end
    
    grain_matching = grain_matching_temp;    
    grain_matching_2 = grain_matching(~isnan(grain_matching(:,2)),:);
    match_eff = size(grain_matching_2,1)/length(numElement_2);
    fprintf(' Assignments made with efficiency %d.\n', match_eff);

end

