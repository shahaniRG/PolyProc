function voi = scope_mask(varargin)

% scope_mask set the volume of interest along different gid_map
%==========================================================================
% FILENAME:          scope_mask.m
% DATE:              1 May, 2019     
% PURPOSE:           Determine common volume of interest
%==========================================================================
%IN :
%    aligned gid_map    : (array) 3D data set of gid_map
%                         no restriction on number of gid_map input
%                         dimension of each gid_map input should match
%
%OUT :
%    voi     : (logical) volume of interest is filled with 1(true)
%                      and the rest of array is filled with 0(false)
%                      The dimension of output is also same as each gid_map
%==========================================================================
%EXAMPLE :
%    gid_map_mask = scope_mask(gid_map_al_t1,gid_map_al_t2);
%==========================================================================

num_data = length(varargin);

voi = varargin{1} ~= 0;

    for i = 2:num_data
        
        voi = voi~=0 & varargin{i}~=0;
        
    end

fprintf('Found intersection volume between datasets.\n');

end