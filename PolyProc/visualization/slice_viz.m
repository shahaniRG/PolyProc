function [] = slice_viz(gid_map, varargin)
% slice_visualization reconstruct specified 2d slice
%==========================================================================
% FILENAME:          slice_visualization.m
% DATE:              1 May, 2020        
% PURPOSE:           reconstruct 2d slice
%==========================================================================
%IN :
%    gid_map   : (array) 3D data set of gid_map
%    z_pos     : slice of interest (along sample's axis of rotation)
%
%OPTIONAL : specify grain coloring criterion 
%     |---no specification  - random coloring (default)
%
%     |---'orientation'     - (string) coloring based on crystallographic dirction of grain
%                             (grain_rodV & crystal structure should be input)
%     |---'size'            - (string) coloring based on grain volume
%                             (numElement should be input)
%     |---'completeness'    - (string) coloring based on completeness of grain
%                             (comp should be input)
%     |---'neighbor'        - (string) coloring based on the number of neighboring grain
%
%     |---'GND'             - (string) coloring based on geometrically necessary dislocation density
%           ('voxel_viz')   - (string) plot extra figure of per pixel based GND density distribution 
%
%OUT :
%    2D reconstruction of a slice
%==========================================================================
%EXAMPLE :
% 0. no specification
%    slice_viz(gid_map_1,'slice',77)
%    slice_viz(gid_map_1,'slice',[70 71 72 73 74 75])
%
% 1. orientation
%    slice_viz(gid_map_1,'slice',77,'orientation', grain_rodV_1, 'cubic')
%    slice_viz(gid_map_1,'slice',[70 71 72 73 74 75],'orientation', grain_rodV_1, 'cubic')
%
% 2. size
%    slice_viz(gid_map_1,'slice',77,'size', numElement_1)
%    slice_viz(gid_map_1,'slice',[70 71 72 73 74 75],'size', numElement_1)
%
% 3. completeness
%    slice_viz(gid_map_1,'slice',77,'completness', comp_1)
%
% 4. neighbor
%    slice_viz(gid_map_1,'slice',77,'neighbor')
%
% 5. Geometrically necessary dislocation(GND) density
%    slice_vis(gid_map_1,'slice',77,'GND', rodV_1)
%
%==========================================================================

%% set up different z_pos
if any(strcmp(varargin,'slice')) 
    idx = find(strcmp(varargin,'slice'))+1;
    z_pos = varargin{idx};
else 
    z_pos = varargin{1};
end

%% set up background

gid_list_tot = unique(gid_map(:));
gid_list_tot = gid_list_tot(2:end);

gid_slice = zeros(size(gid_map,1),size(gid_map,2),length(z_pos));
for j = 1:length(z_pos)
    gid_slice(:,:,j) = gid_map(:,:,z_pos(j));
end

    gid_slice_or = zeros(size(gid_slice,1), size(gid_slice,2), 3);
    gid_slice_tot = zeros(size(gid_slice_or));

gid_list = unique(gid_slice(:));
gid_list = gid_list(2:end);
color_gid_list = zeros(size(gid_list,1), 3);


%% define coloring criteria

if ~any(strcmp(varargin,'orientation') | ...
        strcmp(varargin,'completness') | ...
        strcmp(varargin,'size') | ...
        strcmp(varargin,'neighbor') | ...
        strcmp(varargin,'GND'))
    %random coloring
    color_gid_list = rand(size(gid_list,1),3);
    
elseif any(strcmp(varargin,'orientation')) 
    %coloring based on crystallographic direction
    idx = find(strcmp(varargin,'orientation'))+1;
    grain_rodV = varargin{idx};
    crystal = varargin{idx+1};
    cs = crystalSymmetry(crystal);
    oM = ipfHSVKey(cs);
    oM.inversePoleFigureDirection = zvector;
    for i = 1:length(gid_list)
        q1 = rodrigues2quat(vector3d(grain_rodV((gid_list_tot==gid_list(i)),:)));
        o1 = orientation(q1, cs);
        if isnan(o1.phi1)
            color_gid_list(i,:) = [0 0 0];
        else
            color_gid_list(i,:) = oM.orientation2color(o1);
        end
    end
    
    
elseif any(strcmp(varargin,'completness')) 
    %coloring based on average completeness
    idx = find(strcmp(varargin,'completness'))+1;
    comp = varargin{idx};
    comp_mean = zeros(length(gid_list),1);
    for i = 1:length(gid_list)
        comp_mean(i) = mean(comp(gid_map==gid_list(i)));
    end
    comp_color = jet(101);
    for i = 1:length(gid_list)
        color_gid_list(i,:) = ...
            comp_color(ceil( ( comp_mean(i)-min(comp_mean)+0.001 ) / ( max(comp_mean) -min(comp_mean) ) *100 ),:);
    end
    
    
elseif any(strcmp(varargin,'size')) 
    %coloring based on grain size
    idx = find(strcmp(varargin,'size'))+1;
    numElement = varargin{idx};
    size_color = jet(max(numElement(:,2)));
    for i = 1:length(gid_list)
        color_gid_list(i,:) = size_color(numElement(gid_list_tot==gid_list(i),2),:);
    end

elseif any(strcmp(varargin,'neighbor')) 
    %coloring based on the number of neighbor grain
    adj = imRAG(gid_map);
    adj((adj(:,1) > adj(:,2)),:) = [];
    [~,F] = mode(adj(:));
    neighbor_color = jet(F);
    for i = 1:length(gid_list)
        numNeighbor = sum(sum(adj == gid_list(i)));
        color_gid_list(i,:) = neighbor_color(numNeighbor,:);
    end
    
elseif any(strcmp(varargin,'GND')) 
    idx = find(strcmp(varargin,'GND'))+1;
    rodV = varargin{idx};
%     crystal = varargin{idx+1};
%     cs = crystalSymmetry(crystal);
%     oM = ipfHSVKey(cs);
%     oM.inversePoleFigureDirection = zvector;
    [~,m,n,~]=size(rodV);

    %coloring based on the geometrically necessary dislocation density
    % 1. prompt to define extra inputs - 'cubic', 'fcc/bcc ~~', scale.x scale.y, nu
    % (poisson's ratio - default 0.3)
    % 2. copy paste the GND code
    % 3. sub-section for extra visualization- voxel based (add flag)
    % 4. calculate into per grain base
    % 5. store them a form of color scale
    
    % 1. prompt to define extra inputs
    data_structure_prompts={'Crystal structure (e.g. bcc, fcc, etc)',...
                            'Point Group (e.g. 432)' ,...
                            'Lattice Parameter a (nm)',...
                            'Lattice Parameter b (nm)',...
                            'Lattice Parameter c (nm)',...
                            '\alpha (degrees)',...
                            '\beta (degrees)',...
                            '\gamma (degrees)',...
                            'x axis pixel size (um)',...
                            'y axis pixel size (um)',...
                            'poisson''s ratio' };
    default_data_structure={'bcc',...
                            '432',...
                            '0.2947',...
                            '0.2947',...
                            '0.2947',...
                            '90',...
                            '90',...
                            '90',...
                            '1',...
                            '1',...
                            '0.3'};

        % promt setup
    dims=[1 70];
    opts.Interpreter='tex';
    input_directory=inputdlg(data_structure_prompts,'Define parameters',dims,default_data_structure,opts);
    
        % based on input values, define variables
    cs=crystalSymmetry(input_directory{2},... % point group
                       10*[str2double(input_directory{3}), str2double(input_directory{4}), str2double(input_directory{5})], ... % lattice parameters
                       [str2double(input_directory{6}), str2double(input_directory{7}), str2double(input_directory{8})].*degree); % lattice angles
    dx=str2double(input_directory{9}); % define x-axis scale
    dy=str2double(input_directory{10}); % define y-axis scale

    %a=norm(cs.aAxis); % Lattice Parameter

        % Define Dislocation System and Properties
    if strcmp(input_directory{1},'bcc')
        dS=dislocationSystem.bcc(cs);
    elseif strcmp(input_directory{1},'fcc')
        dS=dislocationSystem.fcc(cs);
    end
        
    nu=str2double(input_directory{11}); 
    dS(dS.isEdge).u=1;
    dS(dS.isScrew).u=1-nu;

    gndslice_tot = zeros(size(gid_slice));
    
    for i = 1:length(z_pos)
        
        slice_ori = orientation(rodrigues2quat(vector3d(squeeze(rodV(1,:,:,z_pos(i))),...
                                                        squeeze(rodV(2,:,:,z_pos(i))),...
                                                        squeeze(rodV(3,:,:,z_pos(i))))));
        
        dSRot=slice_ori*dS;
        
        % 2. calculate GND density
        slice_ori_right = slice_ori(:,[2:end end-1]);
        gX = log(slice_ori_right,slice_ori,'left') ./ dx;
        gX(:,end) = - gX(:,end);
        ori_up = slice_ori([2:end end-1],:);
        gY = log(ori_up,slice_ori,'left') ./ dy;
        gY(end,:) = - gY(end,:);

        %MTEX Retrieved Code (curvature.m)
        % Calculate curvature tensor
        kappa = dyad(gX,tensor([1;0;0])) + ...
                dyad(gY,tensor([0;1;0]));
        kappa{:,3}=NaN;
        kappa = curvatureTensor(kappa,'unit','1/um');

        %Caculate Dislocation Density
        [rho,factor]=fitDislocationSystems(kappa,dSRot);
        alpha=sum(dSRot.tensor .* rho,2);
        alpha.opt.unit='1/um';
        gndvector=factor*sum(abs(rho .* dSRot.u),2);
        gndslice=reshape(gndvector,m,n);
        gndslice_tot(:,:,i) = gndslice;
        
        if any(strcmp(varargin,'voxel_viz')) 
            %visualization
            figure
            imagesc(gndslice)
            c = colorbar;
            CLim=[prctile(gndslice(:),50) prctile(gndslice(:),90)];
            %CLim=[1e12 1e15];
            set(gca,'ColorScale','log');
            set(gca,'CLim',CLim);
            c.Label.String='GND Density (m^{-2})';
            axis off; axis equal; 
        end

        % vizualize in per grain based
        
%         per_grain_slice = zeros(m,n);
%         for k = 1:length(gid_list)
%             temp_slice = gid_map(:,:,z_pos(i));
%             per_grain_slice(temp_slice==(gid_list(k)))= mean(gndslice(temp_slice==(gid_list(k))),'all');            
%         end
%         figure
%         imagesc(per_grain_slice)
%         axis off
    end

%     grain_GND = zeros(size(gid_list,1),2);
%     grain_GND(:,1) = gid_list';
% 
%     for j = 1:length(gid_list)
%         grain_GND(j,2) = mean(gndslice_tot(gid_slice == gid_list(j)));
%     end
%     grain_GND(:,3) = round(log(grain_GND(:,2)).*100);
%     GND_color = jet(max(grain_GND(:,3)));
%     for i = 1:length(gid_list)
%         color_gid_list(i,:) = GND_color(grain_GND(grain_GND(:,1)==gid_list(i),3),:);
%     end   
    

end



%% visualization

if ~any(strcmp(varargin,'GND')) 
for k = 1:length(z_pos)
            gid_slice_layer = gid_slice(:,:,k);
            %change background color to white
            testIm = gid_slice_layer == 0;
            gid_slice_or_copy = gid_slice_or;
            gid_slice_or_copy(:,:,1) =  testIm;
            gid_slice_or_copy(:,:,2) =  testIm;
            gid_slice_or_copy(:,:,3) =  testIm;
            gid_slice_tot_copy = gid_slice_tot + gid_slice_or_copy;


            fprintf('Displaying slice #%d ...\n', z_pos(k));
            gid_list_layer = unique(gid_slice_layer(:));
            gid_list_layer = gid_list_layer(2:end);
            
            for i = 1:length(gid_list_layer)

                testIm = (gid_slice_layer == gid_list_layer(i));
                
                gid_slice_or_copy(:,:,1) = color_gid_list(gid_list==gid_list_layer(i),1) .* testIm;
                gid_slice_or_copy(:,:,2) = color_gid_list(gid_list==gid_list_layer(i),2) .* testIm;
                gid_slice_or_copy(:,:,3) = color_gid_list(gid_list==gid_list_layer(i),3) .* testIm;

                gid_slice_tot_copy = gid_slice_tot_copy + gid_slice_or_copy;
            end

            figure, imagesc(gid_slice_tot_copy);
            axis equal;
            axis off;
            % colorbar setup
                if any(strcmp(varargin,'completness'))
                    colormap(comp_color)
                    caxis([min(comp_mean) max(comp_mean)])
                    colorbar
                elseif any(strcmp(varargin,'size'))
                    colormap(size_color)
                    colorbar
                    caxis([0 max(numElement(:,2))])
                elseif any(strcmp(varargin,'neighbor'))    
                    colormap(neighbor_color)
                    colorbar
                    caxis([0 F]) 
%                 elseif any(strcmp(varargin,'GND'))    
%                     colormap(GND_color)
%                     c = colorbar;
%                     caxis([0 max(grain_GND(:,2))]) 
%                     set(gca,'ColorScale','log');
%                     c.Label.String='GND Density (m^{-2})';
                end
end
end

end
