function Data = cleanup_by_orientation(Data, misoriThresh)
    % Perform misorientation-based clustering and cleanup
    if nargin < 2 || isempty(misoriThresh)
        misoriThresh = 5; % default threshold in degrees
    end

    gid_map = Data.volumeData.uniq_id_phase;
    rodV = Data.volumeData.orientation;
    cs = Data.CrystalSymmetry;

    adj = imRAG(gid_map); adj(adj(:,1) > adj(:,2), :) = [];
    unique_gid = unique(adj);
    label_indices = label2idx(gid_map);
    
    grain_rodV = zeros(max(unique_gid), 3);
    for i = unique_gid(:)'
        idx = label_indices{i};
        grain_rodV(i,:) = [nanmean(rodV(1)(idx)), nanmean(rodV(2)(idx)), nanmean(rodV(3)(idx))];
    end

    q_all = rod2quat(grain_rodV);
    o_all = orientation(q_all, cs);

    adj_misori = false(size(adj,1),1);
    parfor k = 1:size(adj,1)
        adj_misori(k) = angle(o_all(adj(k,1)), o_all(adj(k,2)))/degree < misoriThresh;
    end

    adj_cpy = adj(adj_misori,:);
    A = accumarray(adj_cpy, 1, [max(adj_cpy,[],'all'), max(adj_cpy,[],'all')]);
    C = adj2cluster(A + A');
    
    for i = 1:length(C)
        ids = C{i};
        [~, imax] = max(histc(gid_map(:), ids)); keep_id = ids(imax);
        for j = 1:length(ids)
            gid_map(gid_map == ids(j)) = keep_id;
        end
    end

    Data.volumeData.uniq_id_phase = gid_map;
end
