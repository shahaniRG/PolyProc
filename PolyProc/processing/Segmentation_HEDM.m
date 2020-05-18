function [sample]=Segmentation_HEDM(sample)
% Reconstructs grain from HEDM orientation data imported from Import_HEDM.m
%==========================================================================
% FILENAME:          Segmentation_HEDM.m
% DATE:              19 April, 2020        
% PURPOSE:           Segmentation of grains from orientation data
%==========================================================================
%IN :
%    sample: structure containing orientation, confidence, and scaling data
%            for matrix material and confidence data for particles if 
%            particles are present
%
%OUT :
%    sample.matrix.grain_map: Array of non zero integers assigning
%                             identification numbers to grains
%                                        
%==========================================================================
%EXAMPLE :
%    [Sample1]=Segmentation_HEDM_DEMO(Sample1)
%==========================================================================


%% Establish Dataset Parameters
[m,n,l]=size(sample.matrix.confidence);
No_slices=l;
ori3D=sample.matrix.EulerAngles;
z_space=sample.scale.z;
x_space=sample.scale.x;

%% Get Threshold Values

%Dialgoue Box Formatting
dims=[1 100];
opts.Interpreter='tex';

threshold_prompts={'Grain Boundary Misorientation Threshold (degrees)','Matrix Confidence Threshold'}; 
default_thresholds={'5','0.3'};

thresholds=inputdlg(threshold_prompts,'Thresholds',dims,default_thresholds,opts);

misoriThresh = str2num(thresholds{1});
confiThresh = str2num(thresholds{2});



%% Sample Point Identification
bounds=zeros(m,n,l);
for i=1:No_slices
    
    % Create confidence mask
    confidence=squeeze(sample.matrix.confidence(:,:,i));
    confidence_mask=confidence>confiThresh;
    
    % Get indices of points on the sample
    sample_pts=find(confidence_mask);
    [a,b]=ind2sub([m n],sample_pts);
    ab=[a b];
    
    % Calculate indices of convex hull edge points
    conhull=boundary(a,b,0);
    
    % Get indices of all points
    [e f]=ind2sub([m n],find(ones(m,n)));
    ef=[e f];
    
    % Get convex hull edge points
    conhullpts=[a(conhull),b(conhull)];
    
    % Use inhull to find which points are inside convex hull
    on_sample_pts=inhull(ef,conhullpts);
    on_sample_pts=double(on_sample_pts);
    bounds(:,:,i)=reshape(on_sample_pts,m,n);
end

%% Grain Boundary Identification
gid_mat=zeros(m,n,l);
for i=1:No_slices-1
    % Get sample mask and Euler Angles for slice i
    conhull_mask=squeeze(bounds(:,:,i));
    ori=squeeze(ori3D(:,:,i));
    
    % Calculate misorientation of voxels in x+1, y+1, z+1
    ori_thru=squeeze(ori3D(:,:,i+1));
    ori_right = ori(:,[2:end end-1]);
    ori_up = ori([2:end end-1],:);
    misori_up = angle(ori_up,ori)/degree;
    misori_right=angle(ori_right,ori)/degree;
    misori_thru=angle(ori_thru,ori)/degree;
    
    % Initialize gid map for slice i
    gid_map=conhull_mask;
    
    % Find grain boundaries and set corresponding voxel to -1
    gb=(gid_map~=0) & (misori_up>misoriThresh | misori_right>misoriThresh | misori_thru>misoriThresh);
    gid_map(gb==1)=-1;
    gid_mat(:,:,i)=gid_map;
end

% Calculate grain boundary for last slice

% Get sample mask and Euler Angles for slice i
ori=squeeze(ori3D(:,:,end));
conhull_mask=squeeze(bounds(:,:,end));

% Calculate misorientation of voxels in x+1, y+1
ori_right = ori(:,[2:end end-1]);
ori_up = ori([2:end end-1],:);
misori_up = angle(ori_up,ori)/degree;
misori_right=angle(ori_right,ori)/degree;

% Initialize gid map for slice i
gid_map=conhull_mask;

% Find grain boundaries and set corresponding voxel to -1
gb=(gid_map~=0) & (misori_up>misoriThresh | misori_right>misoriThresh);
gid_map(gb==1)=-1;
gid_mat(:,:,end)=gid_map;

% Create matrix of only grain boundaries
grain_boundaries=double(gid_mat==-1);

% Reassign grain boundaries to 0 for connected component
gid_mat(gid_mat==-1)=0;


%% Assign Grain Identification Numbers

% Use connected components to segment grains
CC=bwconncomp(gid_mat,6);
PixID=CC.PixelIdxList;

% Check each region and possibly assign grain identification number
grain_identification=0;
for j=1:numel(PixID)
    
    % Find indices of region j
    ind=PixID{j};
    
    %Remove region j if volume is less than sphere spanning z_space
        grain_identification=grain_identification+1;
        gid_mat(ind)=grain_identification;
        
    
end

%% Fill Up Grain Boundaries

%Set values of 0 inside sample to NaN
gid_mat(gid_mat==0)=NaN;
gid_mat(bounds==0)=0;

gid_mat=fillmissing(gid_mat,'nearest',2,'EndValues','nearest');
gid_mat=fillmissing(gid_mat,'nearest',1,'EndValues','nearest');

%% Output
sample.matrix.grain_map=gid_mat;
sample.matrix.thresholds.misorientation=misoriThresh;
sample.matrix.thresholds.confidence=confiThresh;

% Remove third dimension for single slice datasets
if No_slices==1
    sample.matrix.EulerAngles=squeeze(sample.matrix.EulerAngles);
    sample.matrix.confidence=squeeze(sample.matrix.confidence);
    sample.matrix.grain_map=squeeze(sample.matrix.grain_map);
end

%% Write to HDF5 File
sampleName=inputname(1);
filename=sprintf('%s.h5',sampleName);
h5create(filename,'/HEDM/Data/GrainId',size(sample.matrix.grain_map));
h5write(filename,'/HEDM/Data/GrainId',sample.matrix.grain_map);
h5create(filename,'/HEDM/Data/Completeness',size(sample.matrix.confidence));
h5write(filename,'/HEDM/Data/Completeness',sample.matrix.confidence);
EulerAngles=zeros(3,m,n,l);
h5create(filename,'/HEDM/Data/EulerAngles',size(EulerAngles));
EulerAngles(1,:,:,:)=ori3D.phi1.*logical(gid_mat);
EulerAngles(2,:,:,:)=ori3D.Phi.*logical(gid_mat);
EulerAngles(3,:,:,:)=ori3D.phi2.*logical(gid_mat);
h5write(filename,'/HEDM/Data/EulerAngles',EulerAngles);
