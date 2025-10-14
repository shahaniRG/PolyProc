function Data = cleanup_by_completeness(Data, completeness_mini)
    % Removes grains whose mean completeness is below threshold
    if nargin < 2 || isempty(completeness_mini)
        completeness_mini = 0.3; % default threshold in degrees
    end
    gid_map = Data.volumeData.uniq_id_phase;
    comp = Data.volumeData.completeness;
    numElement = Data.featureData.numElement;

    label_indices = label2idx(gid_map);
    comp_mean = cellfun(@(idx) mean(comp(idx)), label_indices(1:size(numElement,1)))';
    
    low_comp = numElement(comp_mean < completeness_mini, 1);
    gid_map(ismember(gid_map, low_comp)) = 0;

    Data.volumeData.uniq_id_phase = gid_map;
end
