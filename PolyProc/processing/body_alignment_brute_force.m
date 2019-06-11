function [parameter,gid_map_align,tform,rot_rodri] = body_alignment_brute_force(filename_reference,filename_object,varargin)
 % body_alignment_brute_force automatically align 3D object to register
 % multiple time steps by serial comparison
%==========================================================================
% FILENAME:          body_alignment_brute_force.m
% DATE:              5 Jun, 2019        
% PURPOSE:           3D volume registration by serial comparison
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
%
%OUT :
%    parameter     : (1*6 array) [Rx, Ry, Rz, Tx, Ty, Tz]
%    gid_map_align : (array of 3D dataset) registered grain ID map
%    tform         : (1*1 affine3d) transformation matrix for updating scalar data
%    rot_rodri     : (1*1 vector3d) Euler rotation matrix for updating orientational vector data
%          
%==========================================================================
%EXAMPLE :
%    parameter = body_alignment_brute_force('t1_1.h5','t2_1.h5','angle',[7,7,7]);
%==========================================================================   
tic
    % load data
    Ref = h5read(filename_reference,'/LabDCT/Data/GrainId');
    Target = h5read(filename_object,'/LabDCT/Data/GrainId');
   
    Ref = Ref>0;
    Ref = imfill(Ref,'holes');
    Ref = double(Ref);
    
    Target = Target>0;
    Target = imfill(Target,'holes');
    Target = double(Target);
    
    % defining angular thresholds
    % user can input (1*3) array to specify alpha, beta, gamma, repectively,
    % or input single digit to set alpha, beta, gamma to be same value.
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
        
    % Extend matrix to avoid boundary effect
    extend = round(size(Ref,3)/6);
    RefE = zeros(size(Ref,1)+2*extend,size(Ref,2)+2*extend,size(Ref,3)+2*extend);
    RefE(extend+1:size(Ref,1)+extend,...
        extend+1:size(Ref,2)+extend,extend+1:size(Ref,3)+extend) = Ref;
    Ref = RefE;
    clear RefE 
    
    % Set the range for alignment
    range = round(size(Ref,3)/5);
    fixed = Ref(:,:,round(size(Ref,3)/2)-range:round(size(Ref,3)/2)+range);    
    fixedV = sum(fixed(:));
    fixedArea = sum(fixed,3)>2*range;
    RfixedArea = imref2d(size(fixedArea));
    fixedA = sum(fixedArea(:));
    
    % Run alignment 
    Vvo = [];        
    for a = -angle_x:angle_x        
        for b = -angle_y:angle_y            
            for c = -angle_z:angle_z  
                Rx = [1,0,0,0;0,cosd(a),sind(a),0;0,-sind(a),cosd(a),0;0,0,0,1];
                Ry = [cosd(b),0,-sind(b),0;0,1,0,0;sind(b),0,cosd(b),0;0,0,0,1];
                Rz= [cosd(c),sind(c),0,0;-sind(c),cosd(c),0,0;0,0,1,0;0,0,0,1];
                tform = affine3d(Rx*Ry*Rz);
                Targ = imwarp(Target,tform,'interp','nearest',...
                    'OutputView',imref3d(size(Ref)));
                
                % Find the arropriate fit candidate along z direction
                movV = zeros(1,size(Targ,3)-2*range);
                movA = zeros(1,size(Targ,3)-2*range);
                for i = 1:(size(Targ,3)-2*range)
                    temp = Targ(:,:,i:i+range*2);
                    movV(i) = sum(temp(:));
                    tempArea = sum(temp,3)>2*range;
                    movA(i) = sum(tempArea(:));
                end
                clear i temp tempArea                
                movVR = movV/fixedV;
                movAR = movA/fixedA;
                ind = find(movVR>0.965 & movAR>0.965);
                clear movV movA movVR movAR
                
                % Find the arropriate fit candidate in x-y plane
                if isempty(ind)
                    DR3d = NaN;
                    tx = NaN;
                    ty = NaN;
                    tz = NaN;
                else
                    tempDAR = NaN(1,length(ind));
                    for n = 1:length(ind)
                        indn = ind(n);
                        tempn = Targ(:,:,indn:indn+range*2);
                        tempnArea = sum(tempn,3)>2*range;
                        tformEstimate = imregcorr(tempnArea,fixedArea,'translation');                        
                        tempmovingReg = imwarp(tempnArea,tformEstimate,'OutputView',RfixedArea);
%                         figure
%                         imshowpair(fixed>3,moving>3,'montage')
%                         figure
%                         imshowpair(fixed,movingReg,'montage')                        
                        tempDA = abs(tempmovingReg-fixedArea);
                        tempDAR(n) = sum(tempDA(:))/fixedA;
                    end

                    indm = find(tempDAR==min(tempDAR));
                    indn = ind(indm(1));
                    tempn = Targ(:,:,indn:indn+range*2);
                    tempnArea = sum(tempn,3)>2*range;
                    tformEstimate = imregcorr(tempnArea,fixedArea,'translation'); 
                    tformx = tformEstimate.T(3,1);
                    tformy = tformEstimate.T(3,2);
                    clear tempnArea tformEstimate
                    
                    % Optimize translation in x-y plane
                    Opt = zeros(7,7);
                    for i = 1:7
                        for j = 1:7
                            u = tformx-4+i;
                            v = tformy-4+j;
                            temptformEstimate3d = affine3d([1,0,0,0;0,1,0,0;0,0,1,0;u,v,0,1]);
                            tempmoving3dReg = imwarp(tempn,temptformEstimate3d,'OutputView',...
                                imref3d(size(fixed)));
                            tempDV = abs(tempmoving3dReg-fixed);
                            Opt(i,j) = sum(tempDV(:))/fixedV;                                
                        end
                    end
                    clear i j temptformEstimate3d tempmoving3dReg tempDV u v
                    [x,y] = ind2sub(size(Opt),find(Opt==min(Opt(:))));
                    DR3d = min(Opt(:));
                    tx = x+tformx-4-extend;
                    ty = y+tformy-4-extend;
                    tz = round(size(Ref,3)/2)-indn-extend-range;

                end
                Vvo = [Vvo;a,b,c,DR3d,tx(1),ty(1),tz];
                fprintf('finish %d %d %d....\n',a,b,c)
            end
        end
    end
toc
%     figure
%     plot(Vvo(:,4))
%     figure
%     plot(Vvo(:,5))
%     figure
%     plot(Vvo(:,6))

    goal = find(Vvo(:,4)==min(Vvo(:,4)));
    parameter = zeros(1,6);
    parameter(1:3) = Vvo(goal,1:3);
    parameter(4:6) = Vvo(goal,5:7);
    
    fprintf('%d\n',parameter)
    
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
    
    gid_map_align = imwarp(Target,tform,'OutputView',imref3d(size(Ref)),'interp','nearest');

%% confirmation
    figure
    S_A = isosurface(Ref,0);
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