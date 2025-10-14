function Data = store_absorption_mask_from_workspace(Data)
    % Prompts user to select a variable from the workspace and stores it
    % as Data.volumeData.absorptionMask

    % Prompt user for variable name
    varname = inputdlg('Enter the name of the absorption mask variable in workspace:', ...
                       'Select Absorption Mask', [1 60]);

    if isempty(varname)
        disp('No variable selected. Operation cancelled.');
        return;
    end

    try
        % Evaluate the variable from the base workspace
        mask = evalin('base', varname{1});
    catch
        error('Could not find variable "%s" in workspace.', varname{1});
    end

    % Validate size
    if ~isequal(size(mask), size(Data.volumeData.uniq_id_phase))
        error('Absorption mask must match the size of Data.volumeData.uniq_id_phase.');
    end

    % Ensure logical format
    if ~islogical(mask)
        mask = mask ~= 0;
    end

    % Store into struct
    Data.volumeData.mask = mask;
    disp('Absorption mask successfully stored in Data.volumeData.mask');
end
