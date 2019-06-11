function[FitFun] = reg_solver(a,b,c,x,y,z,Ref,Target)
% This function takes as input three translations and three rotations and
% transforms "next" dataset by these translations; and finally computes the
% fitness (error) between the "prev" and translated "next" datasets.
%==========================================================================
% FILENAME:          reg_solver.m
% DATE:              1 May, 2019     
% PURPOSE:           fitness computation for GA
%==========================================================================
%IN :
%    parameters            : (double) six parameters of transformation
%        alpha, beta, gamma - rotation angle with respect to x-, y-, z-axis
%        x, y, z - translations about x-, y-, z-axis
%    Ref                   : (array of 3D dataset) gid_map of reference
%    Target                : (array of 3D dataset) gid_map of target 
%        to rotate and/or translate to align w.r.t. reference volume
%
%
%OUT :
%    FitFun                : (scalar) % overlap between Ref and Target
%          
%==========================================================================

    rang = round(size(Ref,3)/5);
    prev = Ref(:,:,round(size(Ref,3)/2)-rang:round(size(Ref,3)/2)+rang);
    
% Define the transformation matrix.
    Rx = [1,0,0,0;0,cosd(a),sind(a),0;0,-sind(a),cosd(a),0;0,0,0,1];
    Ry = [cosd(b),0,-sind(b),0;0,1,0,0;sind(b),0,cosd(b),0;0,0,0,1];
    Rz= [cosd(c),sind(c),0,0;-sind(c),cosd(c),0,0;0,0,1,0;0,0,0,1];
    T = [1,0,0,0;0,1,0,0;0,0,1,0;round(x),round(y),round(z),1];
    tform = affine3d(Rx*Ry*Rz*T);
    
% Apply the affine transformation matrix on the "next" dataset
    Targ = imwarp(Target,tform,'interp','nearest',...
              'OutputView',imref3d(size(Ref)));           
    next = Targ(:,:,round(size(Ref,3)/2)-rang:round(size(Ref,3)/2)+rang);

% Compute the fitness function as the fraction of pixels that do not
% overlap between the "prev" and translated next datasets
    temp = prev+next;
    temp = temp>0;
    FitFun1 = sum(temp(:))/sum(prev(:));
    FitFun2 = sum( abs(next(:) - prev(:)) )/sum(prev(:));
    FitFun = FitFun1*FitFun2;

end