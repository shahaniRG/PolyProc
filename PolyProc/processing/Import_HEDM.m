function [sample]=Import_HEDM
% Imports the raw data from a HEDM experiment
%==========================================================================
% FILENAME:          Import_HEDM.m
% DATE:              12 April, 2020        
% PURPOSE:           HEDM Data Importing
%==========================================================================
%IN :
%
%OUT :
%    sample : structure containing orientation, confidence, and scaling data for
%               matrix material and confidence data for particles if
%               particles are present
%                                        
%==========================================================================
%EXAMPLE :
%    [Sample1]=Import_HEDM
%==========================================================================

%% Ask User if Particle are Present
partprompts={'Are secondary phases are present? (1-True, 0-False)'};

%Dialgoue Box Formatting
dims=[1 100];
opts.Interpreter='tex';

part=inputdlg(partprompts,'Particles',dims);
part=str2num(part{1});
%% Get h5 file
matrix_filename=uigetfile('*.h5','Select Matrix Data');

if part==1
particle_filename=uigetfile('*.h5','Select Particle Data');
end

%% Collect User Inputs and Thresholds

%Thresholds

% if part==1
%     threshold_prompts={'Grain Boundary Misorientation Threshold (degrees)','Matrix Confidence Threshold',...
%         'Particle Confidence Threshold'}; 
%     default_thresholds={'5','0.3','0.5'};
% else
%     threshold_prompts={'Grain Boundary Misorientation Threshold (degrees)','Matrix Confidence Threshold'}; 
%     default_thresholds={'5','0.3'};
% end

% Matrix Crystal Symmetry
symmetry_prompts={'Point Group', 'Lattice Parameter a (nm)',...
    'Lattice Parameter b (nm)','Lattice Parameter c (nm)','\alpha (degrees)','\beta (degrees)','\gamma (degrees)'};
default_symmetry={'432','0.2947','0.2947','0.2947', '90','90','90'};

% Particle Crystal Symmetry
if part==1
    partsymmetry_prompts={'Point Group', 'Lattice Parameter a (nm)',...
    'Lattice Parameter b (nm)','Lattice Parameter c (nm)','\alpha (degrees)','\beta (degrees)','\gamma (degrees)'};
    partdefault_symmetry={'432','0.3692','0.3692','0.3692', '90','90','90'};
end

%Scaling
scaling_prompts={'x Spacing (\mum/(pixel or voxel)','y Spacing (\mum/(pixel or voxel)','z Spacing (\mum/(pixel or voxel)'};
default_scaling={'1','1','7'};

%HDF5 Structure
data_structure_prompts={'Directory to Raw Data (Include Slash)','Euler Angle Data Name (Include Slash)',...
    'Confidence Data Name (Include Slash)' };
default_data_structure={'/slices','/EulerAngles','/Confidence'};




matrixsymmetry=inputdlg(symmetry_prompts,'Matrix Crystal Symmetry',dims,default_symmetry,opts);

if part==1
partsymmetry=inputdlg(partsymmetry_prompts,'Particle Crystal Symmetry',dims,partdefault_symmetry,opts);
end

scaling=inputdlg(scaling_prompts,'Scaling',dims,default_scaling,opts);
data_structure=inputdlg(data_structure_prompts,'HDF5 Structure',dims,default_data_structure,opts);

%% Assign Inputs and Thresholds to variables


%Crystal Symmetry
matrix_point_group=matrixsymmetry{1};
matrix_lattice_constants=10*[str2num(matrixsymmetry{2}) str2num(matrixsymmetry{3}) str2num(matrixsymmetry{4})];
matrix_unit_cell_angles=[str2num(matrixsymmetry{5}) str2num(matrixsymmetry{6}) str2num(matrixsymmetry{7})];

if part==1
    part_point_group=partsymmetry{1};
    part_lattice_constants=10*[str2num(partsymmetry{2}) str2num(partsymmetry{3}) str2num(partsymmetry{4})];
    part_unit_cell_angles=[str2num(partsymmetry{5}) str2num(partsymmetry{6}) str2num(partsymmetry{7})];
end

%Scaling
sample.scale.x=str2num(scaling{1});
sample.scale.y=str2num(scaling{2});
sample.scale.z=str2num(scaling{3});


%HDF5 Structure
dir2slices=data_structure{1};
eulerAngleName=data_structure{2};
confidenceName=data_structure{3};

%% Extract Raw Data from HDF5 Structure

%Extract Array size Information
h5_prop=h5info(matrix_filename,dir2slices);
No_slices=numel(h5_prop.Groups);
slice_names={};
for i=1:No_slices
    slice_names{i}=h5_prop.Groups(i).Name;
end
%Check if file names need to be naturally sorted
slice_names=natsortfiles(slice_names);


slice1name=sprintf('%s%s',slice_names{1},confidenceName);
slice1confi=(h5read(matrix_filename,slice1name));

%Array Dimensions
l=No_slices;
[m n]=size(slice1confi);




%% Get Matrix Phase Orientations

%Crystal Symmetry
cs=crystalSymmetry(matrix_point_group,matrix_lattice_constants, matrix_unit_cell_angles*degree);
sample.matrix.cs=cs;
%Initialize
matrixAng=zeros(3,m,n,l);
matrixConfi=zeros(m,n,l);

%Import Angles and Confidence
for i=1:No_slices
    
    %Euler Angles
    sliceAngles=sprintf('%s%s',slice_names{i},eulerAngleName);
    eulerAng = h5read(matrix_filename,sliceAngles);
    ori=orientation.byEuler((reshape(eulerAng(1,:,:),m,n)),(reshape(eulerAng(2,:,:),m,n)),(reshape(eulerAng(3,:,:),m,n)),cs);
    matrixAng(1,:,:,i)=ori.phi1;
    matrixAng(2,:,:,i)=ori.Phi;
    matrixAng(3,:,:,i)=ori.phi2;
    
    %Confidence
    sliceConfidence=sprintf('%s%s',h5_prop.Groups(i).Name,confidenceName);
    conf = h5read(matrix_filename,sliceConfidence);
    matrixConfi(:,:,i)=conf;
end
ori3D=orientation.byEuler(reshape(matrixAng(1,:,:,:),m,n,l),...
                          reshape(matrixAng(2,:,:,:),m,n,l),...
                          reshape(matrixAng(3,:,:,:),m,n,l),cs);

sample.matrix.EulerAngles=ori3D;
sample.matrix.confidence=matrixConfi;



%% Get Particle Phase Orientations
if part==1
    %Initialize
    partcs=crystalSymmetry(part_point_group,part_lattice_constants, part_unit_cell_angles*degree);
    particleAng=zeros(3,m,n,l);
    particleConfi=zeros(m,n,l);

    %Import Angles and Confidence
    for i=1:No_slices
        %Euler Angles
        sliceAngles=sprintf('%s%s',h5_prop.Groups(i).Name,eulerAngleName);
        eulerAng = h5read(particle_filename,sliceAngles);
        ori=orientation.byEuler((reshape(eulerAng(1,:,:),m,n)),(reshape(eulerAng(2,:,:),m,n)),(reshape(eulerAng(3,:,:),m,n)),partcs);
        particleAng(1,:,:,i)=ori.phi1;
        particleAng(2,:,:,i)=ori.Phi;
        particleAng(3,:,:,i)=ori.phi2;
    

        %Confidence
        sliceConfidence=sprintf('%s%s',h5_prop.Groups(i).Name,confidenceName);
        conf = h5read(particle_filename,sliceConfidence);
        particleConfi(i,:,:)=conf;
    end
    ori3Dpart=orientation.byEuler(reshape(particleAng(1,:,:,:),m,n,l),...
                                  reshape(particleAng(2,:,:,:),m,n,l),...
                                  reshape(particleAng(3,:,:,:),m,n,l),partcs);
    sample.particles.confidence=particleConfi;
    sample.particles.EulerAngles=ori3Dpart;
    sample.particles.cs=partcs;
end


