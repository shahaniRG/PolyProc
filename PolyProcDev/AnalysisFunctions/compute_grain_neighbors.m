function Data = compute_grain_neighbors(Data)
    % Computes number of neighbors and neighbor lists for each grain
    % Requires: Data.featureData.id and Data.featureData.adjacency

    grain_ids = Data.featureData.id;
    adj = Data.featureData.adjacency;
    nGrains = numel(grain_ids);

    % Map grain ID → row index
    id_to_idx = containers.Map(grain_ids, 1:nGrains);

    % Initialize
    neighbors_list = cell(nGrains, 1);

    % Fill neighbor lists from adjacency
    for k = 1:size(adj,1)
        id1 = adj(k,1); id2 = adj(k,2);
        if id_to_idx.isKey(id1)
            neighbors_list{id_to_idx(id1)}(end+1) = id2;
        end
        if id_to_idx.isKey(id2)
            neighbors_list{id_to_idx(id2)}(end+1) = id1;
        end
    end

    % Count
    num_neighbors = cellfun(@numel, neighbors_list);

    % Store in struct
    Data.featureData.neighbors_list = neighbors_list;
    Data.featureData.num_neighbors = num_neighbors;
end
