function Data = compute_feature_properties(Data)
    % Computes grain-wise properties and stores them in Data.featureData
    % Fields: id, numElement, adjacency, grain_coord, grain_rodV, comp_avg

    gid_map = Data.volumeData.uniq_id_phase;
    rodV = Data.volumeData.orientation;
    comp = Data.volumeData.completeness;

    % Get list of grain IDs (excluding background ID 0)
    grain_ids = unique(gid_map(:));
    grain_ids(grain_ids == 0) = [];
    nGrains = numel(grain_ids);

    % Initialize feature arrays
    numElement   = zeros(nGrains, 1);
    grain_coord  = zeros(nGrains, 3);
    grain_rodV   = zeros(nGrains, 3);
    comp_avg     = zeros(nGrains, 1);

    label_indices = label2idx(gid_map);
    
    % Build a map from grain ID to label_indices index
    gid_to_idx = containers.Map(grain_ids, 1:nGrains);

    for i = 1:length(label_indices)
        idx = label_indices{i};
        if isempty(idx)
            continue
        end
        gid = gid_map(idx(1));
        if gid == 0 || ~isKey(gid_to_idx, gid)
            continue
        end
        j = gid_to_idx(gid);

        numElement(j) = numel(idx);
        [x, y, z] = ind2sub(size(gid_map), idx);
        grain_coord(j, :) = mean([x, y, z], 1);
        comp_avg(j) = mean(comp(idx));
        grain_rodV(j, :) = [nanmean(rodV(1,idx)), nanmean(rodV(2,idx)), nanmean(rodV(3,idx))];
    end

    % Compute adjacency
    adj = imRAG(gid_map); 
    adj(adj(:,1) > adj(:,2), :) = [];

    % Store in Data.featureData
    Data.featureData.id           = grain_ids;
    Data.featureData.numElement   = numElement;
    Data.featureData.grain_coord  = grain_coord;
    Data.featureData.grain_rodV   = grain_rodV;
    Data.featureData.comp_avg     = comp_avg;
    Data.featureData.adjacency    = adj;
end
