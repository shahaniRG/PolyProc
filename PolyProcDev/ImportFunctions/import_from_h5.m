function Data = import_from_h5(file_name)
% IMPORT_FROM_H5
% Imports volume data (GrainId, Rodrigues, Completeness) from a LabDCT .h5 file
% into the nested struct format: Data.volumeData
% Also prompts the user to select crystal structure

    %% Prompt user to define HDF5 dataset paths
    data_structure_prompts = {
        'Directory for grain ID map', ...
        'Directory for orientation map (Rodrigues)', ...
        'Directory for completeness map'
    };
    default_data_structure = {
        '/LabDCT/Data/GrainId', ...
        '/LabDCT/Data/Rodrigues', ...
        '/LabDCT/Data/Completeness'
    };

    dims = [2 70];
    opts.Interpreter = 'tex';

    data_directory = inputdlg(data_structure_prompts, 'HDF5 Structure', dims, default_data_structure, opts);

    %% Read data from file into struct
    Data.volumeData.uniq_id_phase  = h5read(file_name, data_directory{1});
    Data.volumeData.orientation    = h5read(file_name, data_directory{2});
    Data.volumeData.completeness   = double(h5read(file_name, data_directory{3}));

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
