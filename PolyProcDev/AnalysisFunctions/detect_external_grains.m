function Data = detect_external_grains(Data)
    % Marks grains that touch air as surface grains
    % Requires: Data.volumeData.uniq_id_phase, Data.volumeData.mask
    % Stores:   Data.featureData.is_surface (logical array)

    gid_map  = Data.volumeData.uniq_id_phase;
    raw_mask = Data.volumeData.mask;
    grain_ids = Data.featureData.id;

    if isempty(raw_mask)
        error('Mask not found in Data.volumeData.mask.');
    end

    % --- Fill holes in mask ---
    filled_mask = imfill(raw_mask, 'holes');

    % --- Find surface voxels (mask edge) ---
    se = strel('sphere', 1);
    mask_eroded = imerode(filled_mask, se);
    boundary_voxels = filled_mask & ~mask_eroded;

    % --- Grain IDs touching the boundary ---
    external_voxels = gid_map(boundary_voxels);
    surface_ids = unique(external_voxels);
    surface_ids(surface_ids == 0) = [];  % exclude background

    % --- Create flag per grain ID ---
    id_to_surface = ismember(grain_ids, surface_ids);
    Data.featureData.is_surface = id_to_surface;

    fprintf('[detect_external_grains] Found %d surface grains.\n', sum(id_to_surface));
end
