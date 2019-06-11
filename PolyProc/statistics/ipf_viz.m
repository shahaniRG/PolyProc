function [] = ipf_viz(gid_map, grain_rodV, crystal)
% ipf_visualization generates ipf with directions of every grains
%==========================================================================
% FILENAME:          ipf_visualization.m
% DATE:              1 May, 2019     
% PURPOSE:           generate IPF
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%    grain_rodV: n*3 array with rodrigues vectors of grains
%    crystal   : (string) 'cubic' for cubic crystal structure
%
%OUT :
%    Inverse Pole Figure
%       ( grain color is based on normal direction of grain)
%==========================================================================
%EXAMPLE :
%           ipf_viz(gid_map_2, grain_rodV_2, cs_2)
%==========================================================================
cs = crystalSymmetry(crystal);
gid_list = unique(gid_map(:));
gid_list = gid_list(2:end);
oM = ipfHSVKey(cs);
oM.inversePoleFigureDirection = zvector;

figure;
for i = 1:length(gid_list)
    
    q1 = rodrigues2quat(vector3d(grain_rodV(i,:)));
    o1 = orientation(q1, cs);
    
    plotIPDF(o1,oM.orientation2color(o1),zvector,'MarkerSize',10);
    hold on;
    
end

end