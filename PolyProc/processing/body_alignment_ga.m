function [parameter,gid_map_align,tform,rot_rodri] = body_alignment_ga(filename_reference,filename_object,varargin)
% body_alignment_auto automatically align 3D object to register multiple time steps
%==========================================================================
% FILENAME:          body_alignment_ga.m
% DATE:              1 May, 2020      
% PURPOSE:           automated 3D volume registration
%==========================================================================
%IN :
%    filename_reference    : (string) name of hdf5 file that you are trying
%                               to use as reference
%    filename_object       : (string) name of hdf5 file taht you are trying
%                               to rotate and/or translate to align with reference object
%
%OPTIONAL :
%    angle      : setting up the range of rotational angle.
%                 (double) This value applies to all alpha, beta, gamma roation angle range
%                   default : 5
%                 (1*3 array) each different range of alpha, beta, gamma rotation angle are available
%    threshold  : (double) setting up the misfit tolerance thershold
%                   defulat : 0.02
%    visualization : visualize two volumes to confirm alignment
%                   (string) string flag 'viz' will trigger visualization
%                   default : off
%
%OUT :
%    parameter     : (1*6 array) [alpha, beta, gamma, a, b, c]
%    gid_map_align : (array of 3D dataset) registered grain ID map
%    tform         : (1*1 affine3d) transformation matrix for updating scalar data
%    rot_rodri     : (1*1 vector3d) Euler rotation matrix for updating orientational vector data
%          
%==========================================================================
%EXAMPLE :
%    Parameter = body_alignment_ga('t1_1.h5','t2_1.h5','angle',[7,7,7],'viz');
%    (out) Parameter = [0.44 , 6.41 , -1.11, 5.59, 4.44, 15.60]
%==========================================================================
tic
    Ref = h5read(filename_reference,'/LabDCT/Data/GrainId');
    Target = h5read(filename_object,'/LabDCT/Data/GrainId');

    Ref_ori = Ref;
    Target_ori = Target;
    
    Ref = Ref>0;
    Ref = imfill(Ref,'holes');
    Ref = double(Ref);

    Target = Target>0;
    Target = imfill(Target,'holes');
    Target = double(Target);
    
% defining angular thresholds
% user can input (1*3) array to specify alpha, beta, gamma, repectively,
%       or input single digit to set alpha, beta, gamma to be same value.
        if	any(strcmp(varargin,'angle'))
            idx = find(strcmp(varargin,'angle'))+1;
            angle= varargin{idx};
            if length(angle)<3
                angle_x = angle;
                angle_y = angle;
                angle_z = angle;
            else
                angle_x = angle(1);
                angle_y = angle(2);
                angle_z = angle(3);
            end
        else
            angle_x = 5;
            angle_y = 5;
            angle_z = 5;
        end
% Constrain the available translations to the image domain.  That is, only
% allow those translations that do not move the sample outside of the
% image domain, defined by its size in pixels.

    rangx = size(Target,1)/2;
    rangy = size(Target,2)/2;
    rangz = size(Target,3)/2;
    
    ub = [angle_x,angle_y,angle_z,rangx,rangy,rangz];
    lb = [-angle_x,-angle_y,-angle_z,-rangx,-rangy,-rangz];

% Specify some parameters needed to run GA.  Depending on the problem, some
% of these parameters may require some tweaking.  For details on each, see 
% https://www.mathworks.com/help/gads/genetic-algorithm-options.html

% specify fitness threshold based on volex matching
        if	any(strcmp(varargin,'threshold'))
            idx = find(strcmp(varargin,'threshold'))+1;
            threshold= varargin{idx};
        else
            threshold = 0.02;
        end

    fittol = threshold;        % fitness tolerance
    rat = sum(Target(:))/sum(Ref(:));
    numsample = round(rat)*150;
    opts = optimoptions('ga','PopulationSize', numsample, ...
        'FitnessLimit', fittol, 'MaxGenerations',25, 'PlotFcn',@gaplotbestf);

% Run the GA with a simultaneous visualization of the result.  
    [x,fitfun] = ga(@(x)reg_solver(x(1),x(2),x(3),x(4),x(5),x(6), ...
        Ref,Target),6,[],[],[],[],lb,ub,[],opts);
    bestx(1,:) = x;       % Set the best translations to be the ones just found
    beste(1) = fitfun;    % Set the best fitness to be the one just found
    parameter = bestx;
% The "best" translations just determined may be located inside a local 
% minimum that may be still above our required tolerance level. So, keep 
% running the GA until you find a set of translations that is below fittol.

    while fitfun > fittol
    
        fprintf('Misfit tolerance is not satisfied in first iteration.\n Initiating second iteration.\n')
    % Run the GA with a simultaneous visualization of the result.  
        opts = optimoptions('ga','PopulationSize', numsample*2, ...
                  'FitnessLimit', fittol, 'MaxGenerations',25, 'PlotFcn',@gaplotbestf);

        [x,fitfun] = ga(@(x)reg_solver(x(1),x(2),x(3),x(4),x(5),x(6), ...
                          Ref,Target),6,[],[],[],[],lb,ub,[],opts);
    
    % If the new translations are better than the last, store in array
        if fitfun < beste 
            bestx = x;
            beste = fitfun;     
            parameter = x;   
            break        
        end
    
        bestx(size(bestx,1)+1,:) = x; 
        beste(size(beste)+1) = fitfun;        
        ind = beste==min(beste);
        parameter = bestx(ind,:);    
        break
    end
    
    if min(beste) < fittol
        fprintf('Misfit tolerance is satisfied.\n Ouptting parameters.\n')
    elseif min(beste) >= fittol
        fprintf('Misfit tolerance is not satisfied after two iterations.\n Ouptting best-fitting parameters.\n')
    end
toc    

%% Generate transformation and Roatation function using result

Rx = parameter(1);
Ry = parameter(2);
Rz = parameter(3);
Tx = parameter(4);
Ty = parameter(5);
Tz = parameter(6);

% Rotation
    % Rotation about x axis
    Rotation_x = [1,0,0,0;0,cosd(Rx),sind(Rx),0;0,-sind(Rx),cosd(Rx),0;0,0,0,1];
    % Rotation about y axis
    Rotation_y = [cosd(Ry),0,-sind(Ry),0;0,1,0,0;sind(Ry),0,cosd(Ry),0;0,0,0,1];
    % Rotation about y axis
    Rotation_z= [cosd(Rz),sind(Rz),0,0;-sind(Rz),cosd(Rz),0,0;0,0,1,0;0,0,0,1];
    
% Translation
    Translation = [1,0,0,0;0,1,0,0;0,0,1,0;Tx,Ty,Tz,1];
    
% Euler rotation matrix

    rot_x = rotation('axis', xvector,'angle', -Ry*degree);
    rot_y = rotation('axis', yvector,'angle', -Rx*degree); % isosurface flips x & y axis
    rot_z = rotation('axis', zvector,'angle', -Rz*degree); % %mtex and matlab axis direction is opposite
    rot_euler = rot_x*rot_y*rot_z;
    rot_rodri = Rodrigues(rot_euler);

% Full transformation matrix
    
    tform = affine3d(Rotation_x*Rotation_y*Rotation_z*Translation);
    
% apply the alignment
    
    gid_map_align = imwarp(Target_ori,tform,'OutputView',imref3d(size(Ref)),'interp','nearest');
    
%% confirmation
    if	any(strcmp(varargin,'viz'))
    
        figure
        S_A = isosurface(Ref_ori,0);
        patch(smoothpatch(S_A,1,3), 'Facecolor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.1);    
        hold on
        S_B = isosurface(gid_map_align,0);
        patch(smoothpatch(S_B,1,3), 'Facecolor', 'red', 'EdgeColor', 'none','FaceAlpha',0.3);    
        hold off
        xlabel('x');
        ylabel('y');
        zlabel('z');
        title_name = sprintf('%s %s alignment',filename_object,filename_reference);
        title(title_name)
        axis equal
        view([20,30])
        legend(filename_reference,filename_object)
    end

end
