    function [numNeighbor, mis_area, grain_info] = grain_stat(gid_map, numElement,adj,grain_rodV,grain_surface,crystal)
% grain_stat generate figures of four statistical results: histograms of 
% grain volume, topology (number of neighbors), misorientation of grain 
% boundaries, and sphericity. Currently only supports cubic symmetries.
%==========================================================================
% FILENAME:          grain_stat.m
% DATE:              1 May, 2020        
% PURPOSE:           statistical results of grains and their boundaries
%==========================================================================
%IN :
%    gid_map    : (array) 3D data set of gid_map
%
%    numElement : n*2 array with (1st column-gid) & (2nd column-number of voxels)
%
%    adj        : n*2 array shows grain adjacency in number ascending order
%    
%    grain_rodV : n*3 array with rodrigues vectors of grains
%
%    grain_surface : n*1 array wigh grain ID of surface touching grains
%
%    crystal  : (string) crystal structure
%
%OUT :
%    numNeighbor     : n*2 array shows the number of neighboring grains of each grain
%
%    mis_area        : n*4 array with first two column is identical to adj array.
%                      3rd column is misorientation between two adjcent grain.
%                      4th columen is area of grain boundary by two neighboring grains.
%  
%    grain_info      : n*5 array with columns: grain id, number of voxels
%                      per grain, grain surface area, grain volume, grain
%                      sphericity.  Each row is a different grain.
%==========================================================================
%EXAMPLE :
%    [numNeighbor, mis_area, grain_info] = grain_stat(gid_map_1, numElement_1,adj_1,grain_rodV_1,grain_surface,'cubic');
%==========================================================================


%% Distribution of Grain Size

fprintf('Plotting distribution of grain size.\n');

equi_radius = (.75/pi()*numElement(:,2)).^(1/3);
figure(1);
h5 = histogram(equi_radius);
h5.Normalization = 'pdf';
title('Grain size statistics')
xlabel('equivalent grain radious (\mum)')
ylabel('probability density \itp(x)')

hold on
pd = fitdist(equi_radius,'lognormal');
x_values = min(equi_radius)-2:(max(equi_radius)-min(equi_radius))/100:max(equi_radius)+2;
y = pdf(pd,x_values);
plot(x_values,y,'--','LineWidth',2)
legend('experimental data','log-normal distribution')

%% Distribution of Neighbors

fprintf('Plotting distribution of grain neighbors.\n');
numNeighbor = numElement;
numNeighbor_interior = numElement;
numNeighbor_exterior = numElement;

for i = 1:length(numElement)
    
    numNeighbor(i,2) = sum(sum(adj==numElement(i,1)));
    %differentiate exterior grain vs. interir grain
    if ~ismember(numElement(i,1),grain_surface)
        numNeighbor_interior(i,2) = sum(sum(adj==numElement(i,1)));
        numNeighbor_exterior(i,2) = nan;
    else
        numNeighbor_interior(i,2) = nan;
        numNeighbor_exterior(i,2) = sum(sum(adj==numElement(i,1)));
    end
end

h_int = histcounts(numNeighbor_interior(:,2),'BinWidth',2);
h_ext = histcounts(numNeighbor_exterior(:,2),'BinWidth',2);

h_tot = zeros( max([length(h_int) length(h_ext)]) , 2);
h_tot(1:length(h_int),1) = h_int';
h_tot(1:length(h_ext),2) = h_ext';

figure(2)
bar(h_tot, 'stacked')
xlabel('Number of neighbors')
ylabel('Frequency (number of grains)')
title('Topology statistics')
legend('interior grains', 'exterior grains');

%% Misorientation distribution based on Grain Boundary Voxels

fprintf('Plotting distribution of grain boundary misorientation.\n');
mis_area = adj;
cs = crystalSymmetry(crystal);

%calculating (misorientation & area) of each grain boundary
for j = 1:length(adj)
    
    fprintf(' Calculating misorientation of grain boundary #%d...\n',j);
    
    grain_1_idx = numElement(:,1)==adj(j,1);
    grain_2_idx = numElement(:,1)==adj(j,2);
    q1 = rodrigues2quat(vector3d(grain_rodV(grain_1_idx,:)));
    o1 = orientation(q1, cs);
    q2 = rodrigues2quat(vector3d(grain_rodV(grain_2_idx,:)));
    o2 = orientation(q2, cs);
    misori = angle(o1,o2)/degree;
    mis_area(j,3) = misori;
    
    gidK = gid_map ~= adj(j,1) & gid_map ~= adj(j,2);
    gidLocal = double(gid_map);
    gidLocal(gidK) = NaN;
    
    FV = isosurface(gidLocal, (adj(j,1)+adj(j,2))/2);
    
    if size(FV.vertices,1) > 5
        
        FV = smoothpatch(FV,1,3);
        verts = FV.vertices;
        faces = FV.faces;
        e21 =  verts(faces(:, 2), :) - verts(faces(:, 1), :);
        e31 = verts(faces(:, 3), :) - verts(faces(:, 1), :);
        each_area = cross(e21, e31, 2);
        mis_area(j,4) = 1/2 * sum(sqrt(sum(each_area.^2, 2)));
        
    else 
        mis_area(j,4) = 0;
    end
end

figure(3)
mis_area_plot = zeros(round(max(mis_area(:,3))),2);
for k = 1:(round(max(mis_area(:,3)))-1)
    idx = mis_area(:,3)>k & mis_area(:,3)<k+1;
    mis_area_plot(k,1) = k+0.5;                 %first colume is interval of misorientation
    mis_area_plot(k,2) = sum(mis_area(idx,4));  %second colume is sum of boundary area
end    
mis_area_plot(:,3) = mis_area_plot(:,2)./sum(mis_area_plot(:,2)); %third colume is fraction of area

% set up binning to have same interval as background Mackanzie plot 
x_axis = linspace(0, 70, 20);
[h, ~, ~, ~] = histcn(mis_area_plot(:,1), x_axis, ...
        'AccumData',mis_area_plot(:,3));
bar(x_axis(1:end-1),h','barwidth',0.3)

    %With Mackenzie Distribution
        cs = crystalSymmetry('cubic');
        %Inputs of Binwidth and Number of Random numbers
            rndnum = 5000;
        %Create Random Euler Matrix    
            orient_rand = zeros(rndnum,1);
            orient_rand(:,1) = 2*pi()*rand(rndnum,1);
            orient_rand(:,2) = pi()*rand(rndnum,1);
            orient_rand(:,3) = 2*pi()*rand(rndnum,1);
        %Pre-Allocate Matrix and k
            misori_2 = zeros(rndnum,1);
            area = ones(rndnum,1);
        %Apply Crystal Symmetry to find rotation angle
            for i = 1:rndnum
                euler = orient_rand(i,:);
                o1 = orientation('Euler',euler, cs);
                o2 = [0,0,0];
                misori_2(i) = angle(o1,o2)/degree;
            end
    hold on   
    [h2, ~, ~, ~] = histcn(misori_2, x_axis, ...
        'AccumData',area./sum(area(:)));
        bar(x_axis(1:end-1),h2','FaceAlpha',0.2);

title('Misorientation distribution')
ylabel('Grain boundary area fraction')
xlabel('Misorientation (Deg)')
legend('Input data','Random')
axis tight

%% Sphericity
% fprintf('Plotting distribution of grain sphericity.\n');
% for i = 1:length(numElement)
%     
%     FV = isosurface(gid_map == numElement(i,1));
%     FV = smoothpatch(FV,1,3);
%     verts = FV.vertices;
%     faces = FV.faces;
%     if size(FV.vertices,1) > 0
%         e21 =  verts(faces(:, 2), :) - verts(faces(:, 1), :);
%         e31 = verts(faces(:, 3), :) - verts(faces(:, 1), :);
%         each_area = cross(e21, e31, 2);
%         area = 1/2 * sum(sqrt(sum(each_area.^2, 2)));
% 
%         volume = sum(sum(sum(gid_map == numElement(i,1))));
% 
%         numElement(i,3) = area;
%         numElement(i,4) = volume;
%     end
% end
% 
% numElement(:,5) = (pi()^(1/3).*(6*numElement(:,4)).^(2/3))./(numElement(:,3));
% grain_info = numElement;
% 
% figure(4)
% histogram(numElement(:,5),20,'EdgeColor','k');
% xlabel('Sphericity')
% ylabel('Frequency (number of grains)')
% title('Sphericity statistics')
end
