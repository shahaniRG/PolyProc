function Data = import_from_workspace()
% IMPORT_FROM_WORKSPACE
% Imports selected volume data (Grain ID Map, 4D Rodrigues vectors, Completeness Map) 
% from variables currently in the MATLAB workspace into a nested struct format: Data.volumeData

    %% Field display names for user prompt
    field_display = {
        'Grain ID Map', ...
        '4D Rodrigues vectors', ...
        'Completeness Map'
    };
    
    field_keys = {
        'uniq_id_phase', ...
        'orientation', ...
        'completeness'
    };

    imported_fields = struct();

    %% Loop through each field and prompt user for import
    for i = 1:numel(field_keys)
        field_prompt = sprintf('Would you like to import the %s from the workspace?', field_display{i});
        choice = questdlg(field_prompt, ...
                          'Import from Workspace', ...
                          'Yes', 'No', 'Yes');
        if strcmp(choice, 'Yes')
            % Use listdlg to select variable from base workspace
            vars = evalin('base', 'who');
            [indx, tf] = listdlg('PromptString', ...
                                 sprintf('Select the variable corresponding to: %s', field_display{i}), ...
                                 'SelectionMode', 'single', ...
                                 'ListString', vars, ...
                                 'ListSize', [300, 400]);  % Larger dialog box size
            if tf
                var_name = vars{indx};
                imported_fields.(field_keys{i}) = evalin('base', var_name);
            else
                warning('%s not selected. Skipping.', field_display{i});
            end
        end
    end

    %% Assign selected fields to Data.volumeData
    Data.volumeData = imported_fields;

    %% Prompt user for crystal structure
    crystal_options = {'cubic', 'hexagonal', 'orthorhombic', 'tetragonal', 'monoclinic', 'triclinic'};
    [indx, tf] = listdlg('PromptString', 'Select the crystal structure:', ...
                         'SelectionMode', 'single', ...
                         'ListString', crystal_options);

    if tf
        Data.CrystalStructure = crystalSymmetry(crystal_options{indx});
    else
        error('Crystal structure selection cancelled.');
    end
end
