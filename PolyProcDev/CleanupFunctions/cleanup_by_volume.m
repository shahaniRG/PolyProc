function Data = cleanup_by_volume(Data, grain_size_mini)
    % Removes grains below volume threshold
    if nargin < 2 || isempty(grain_size_mini)
        grain_size_mini = 10; % default threshold in degrees
    end

    gid_map = Data.volumeData.uniq_id_phase;
    valid = gid_map(:); valid(valid==0) = [];
    numElement = accumarray(valid,1);
    
    too_small = find(numElement < grain_size_mini);
    gid_map(ismember(gid_map, too_small)) = 0;
    
    % Keep valid elements
    numElement(too_small) = [];
    valid_ids = find(numElement);
    
    Data.volumeData.uniq_id_phase = gid_map;
end
