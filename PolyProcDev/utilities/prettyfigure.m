function [] = prettyfigure(variable,steps)
% inherited from Paul Chao summer 2022
% v0.2(10/12/2022)
% v0.3(12/4/2023) added option to increase smoothing iterations rather than
% rotating. uses same gif scripts
%%
fps = 60;
dt = 1/fps;

%must have figure opened as 'Figure 1' in matlab
Fig = gcf;

set(gcf,'color','w');
axis vis3d
axis off
grid off
% variable = 'figureB';
filename1 = [variable '_1.gif'];
% filename2 = [variable '_1_alphashape3.gif'];
% filename3 = [variable '_1_convexhull3.gif'];

rotation = linspace(0,359,steps);
%rotation = linspace(0,360,2);
n=1;
 for i = rotation
     delete(findall(Fig,'Type','light'))

     campos([0,10,0]);
     camtarget([0,10,10]);

     view([i 18]);hold on,
    delete(findall(gcf,'Type','light')), camlight right, material dull, hold off
     % camlight right
     % camlight left

     %direction = [1 1 0];
     %rotate(plt,direction,10)
     drawnow
     %pause(0.5)
     
      % Capture the plot as an image 
      frame1 = getframe(Fig); 
      im1 = frame2im(frame1); 
      [imind1,cm1] = rgb2ind(im1,256); 
      
%       frame2 = getframe(Fig2); 
%       im2 = frame2im(frame2); 
%       [imind2,cm2] = rgb2ind(im2,256); 
%       
%       frame3 = getframe(Fig3); 
%       im3 = frame2im(frame3); 
%       [imind3,cm3] = rgb2ind(im3,256); 
    %imwrite(im1,[variable '_1_scatter3_' num2str(i)  '.tif']); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind1,cm1,filename1,'gif', 'Loopcount',inf,'DelayTime',dt); 
          %imwrite(imind2,cm2,filename2,'gif', 'Loopcount',inf); 
          %imwrite(imind3,cm3,filename3,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind1,cm1,filename1,'gif','WriteMode','append','DelayTime',dt); 
         % imwrite(imind2,cm2,filename2,'gif','WriteMode','append'); 
         % imwrite(imind3,cm3,filename3,'gif','WriteMode','append'); 
      end 
      n = n+1;
    % saveas(im1,[variable '_1_scatter3_' num2str(i)  '.tif'])
     %saveas(Fig1,[variable '_2_alphashap_' num2str(i)  '.tif'])
     %saveas(Fig2,[variable '_3_convexhull_' num2str(i)  '.tif'])

     % fprintf(strcat('writing frame for rotation: ',num2str(i),'\n'))
     
 end