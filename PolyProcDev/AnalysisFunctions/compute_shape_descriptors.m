function Data = compute_shape_descriptors(Data)
    % Computes geometric shape descriptors for each grain
    % Includes: convexity, circularity, sphericity, elongation, anisotropy,
    %           bounding box ratio, PCA axis lengths

    gid_map   = Data.volumeData.uniq_id_phase;
    grain_ids = Data.featureData.id;
    volumes   = Data.featureData.numElement;
    nGrains   = numel(grain_ids);

    % Preallocate
    convexity         = zeros(nGrains, 1);
    circularity       = zeros(nGrains, 1);
    sphericity        = zeros(nGrains, 1);
    elongation        = zeros(nGrains, 1);
    anisotropy_index  = zeros(nGrains, 1);
    bounding_box_ratio = zeros(nGrains, 1);
    axis_lengths      = zeros(nGrains, 3);  % major, intermediate, minor

    bw = false(size(gid_map));

    for i = 1:nGrains
        id = grain_ids(i);
        bw(:) = gid_map == id;

        % Extract voxel coords
        [x, y, z] = ind2sub(size(gid_map), find(bw));
        coords = [x, y, z];

        % Skip empty
        if isempty(coords)
            convexity(i) = NaN; circularity(i) = NaN; sphericity(i) = NaN;
            axis_lengths(i,:) = NaN; elongation(i) = NaN;
            anisotropy_index(i) = NaN; bounding_box_ratio(i) = NaN;
            continue;
        end

        % Surface
        [faces, verts] = isosurface(bw, 0.5);
        if isempty(faces)
            convexity(i) = NaN; circularity(i) = NaN; sphericity(i) = NaN;
        else
            % Surface area
            face_area = triangle_areas(verts, faces);
            S = sum(face_area(:));

            % Convex hull
            try
                [~, cv_vol] = convhull(verts);
                convexity(i) = volumes(i) / cv_vol;
            catch
                convexity(i) = NaN;
            end

            % Compactness (circularity-like)
            circularity(i) = (volumes(i)^2) / (S^3 + eps);

            % Sphericity
            sphericity(i) = (pi^(1/3) * (6*volumes(i))^(2/3)) / (S + eps);
        end

        % PCA axis lengths
        coords = coords - mean(coords, 1);
        C = cov(coords);
        eigvals = sort(eig(C), 'descend');
        axis_lengths(i,:) = sqrt(eigvals);

        % Elongation & Anisotropy
        elongation(i) = 1 - axis_lengths(i,3) / axis_lengths(i,1);
        anisotropy_index(i) = axis_lengths(i,1) / axis_lengths(i,3);

        % Bounding box volume
        bbox_vol = range(x)*range(y)*range(z);
        bounding_box_ratio(i) = volumes(i) / (bbox_vol + eps);
    end

    % Store everything
    SF = Data.featureData.shape_factors;
    SF.convexity         = convexity;
    SF.circularity       = circularity;
    SF.sphericity        = sphericity;
    SF.axis_lengths      = axis_lengths;
    SF.elongation        = elongation;
    SF.anisotropy_index  = anisotropy_index;
    SF.bounding_box_ratio = bounding_box_ratio;
    Data.featureData.shape_factors = SF;
end
