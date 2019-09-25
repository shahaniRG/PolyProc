function [gid_map]=slice_clean_up(filename,slice_position,cs,misoriThresh,confiThresh)
% Identifies, segments, and maps grains for 2D slices of HEDM data
%==========================================================================
% FILENAME:          slice_clean_up.m
% DATE:              25 September, 2019      
% PURPOSE:           Grain segmentation
%==========================================================================
%IN :
%   filename :       (string) name of hdf5 file of interest
%                               
%   slice_position : (string) location of the slice of interest
%                               
%   cs :             (string) crystal symmetry (default: 'cubic')   
%
%   misoriThresh :   (number) misorientation threshold angle for identifying grain
%                     boundaries (degrees) (default: 4)
%
%   confiThresh :    (number) confidence threshold with a range of [0 1] used
%                     to identify off sample vs. on sample points (default: 0.3)
%OPTIONAL :
%   cs :             (string) crystal symmetry (default: 'cubic')   
%
%   misoriThresh :   (number) misorientation threshold angle for identifying grain
%                     boundaries (degrees) (default: 4)
%
%   confiThresh :    (number) confidence threshold with a range of [0 1] used
%                     to identify off sample vs. on sample points (default: 0.3)
%  
%OUT :
%   gid_map :        (double) matrix corresponding pixel index to grain ID
%          
%==========================================================================
%EXAMPLES :
%    In :  [gid_map]=slice_clean_up('sample1_rt_900_furnace_nf_copperBCC_q11_rot720_37layers_1050x1050_0.001_shift_0_0_0.h5','z19','432',4,0.3);
%    Out : (figure window with color coded grains)

%    In :  [gid_map]=slice_clean_up('sample1_rt_900_furnace_nf_copperBCC_q11_rot720_37layers_1050x1050_0.001_shift_0_0_0.h5','z19','432',4,0.3)
%    Out : (figure window with color coded grains and matrix gid_map)
%==========================================================================

if nargin==2
    cs='cubic'
    misoriThresh=4
    confiThresh=0.4
end

% Import Data
eularAng = h5read(filename,sprintf('/slices/%s/EulerAngles',slice_position));
confidence = h5read(filename,sprintf('/slices/%s/Confidence',slice_position));
gid_map = reshape(1:1:(size(eularAng,2)*size(eularAng,3)), size(eularAng,2),size(eularAng,3));
[l,m,n]=size(eularAng);
ori=orientation.byEuler((reshape(eularAng(1,:,:),m,n)),(reshape(eularAng(2,:,:),m,n)),(reshape(eularAng(3,:,:),m,n)),cs);

% Confidence Mask
confidence_mask = confidence > confiThresh;
gid_map = gid_map.*double(confidence_mask);

% Compute Misorientation
ori_right = ori(:,[2:end end-1]);
ori_up = ori([2:end end-1],:);
misori_up = angle(ori_up,ori)/degree;
misori_right=angle(ori_right,ori)/degree;

% Identify Grain Boundaries
gb=(gid_map~=0) & (misori_up>misoriThresh | misori_right>misoriThresh);
gid_map(gb==1)=0;

% Segment Grains
CC=bwconncomp(gid_map,4);
PixID=CC.PixelIdxList;
No_grains=numel(PixID);

% Define Orientation Color Key
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = zvector;
color_gid_list=oM.orientation2color(ori);
color_av=zeros(No_grains,1);
gid_map_color=zeros(m,n);

% Index Grain IDs and Average Color Values
for i=1:No_grains
    ind=PixID{i};
    gid_map(ind)=i;
    color_av(i)=mean(color_gid_list(ind));
    gid_map_color(ind)=color_av(i);
end

% Visualization
gid_map_color=flip(gid_map_color);
im=imagesc(gid_map_color);
xlim([1 m]);
ylim([1 n]);
axis square
