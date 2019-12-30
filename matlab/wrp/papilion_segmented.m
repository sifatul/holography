addpath('D:\Mathlab\wrp\data') 
close all;
obj_name = 'papilon';

pc = pcread (strcat(obj_name,'.ply'));
pcshow(pc);xlabel('X');ylabel('Y');zlabel('Z');

roi = [  1 512  1   512  205 255 ]; %% 165 
pc_top = select(pc,findPointsInROI(pc,roi));
figure; pcshow(pc_top)

obj_top = im2double(pc_top.Location);
obj_top(:,3) = obj_top(:,3)+ 50;
% obj_c_top = im2double(pc_top.Color); 

pc_new = pointCloud(obj_top);
pc_new.Color=pc_top.Color
figure,pcshow(pc_new); xlabel('X'); ylabel('Y');zlabel('Z');



%% middle

roi = [  1 512  1   512  100 204 ]; %% 165 
pc_middle = select(pc,findPointsInROI(pc,roi));
figure;pcshow(pc_middle); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% obj_last = im2double(pc_last.Location);
% obj_last(:,3) = min(obj_last(:,3));
% 
% 
% pc_new_last = pointCloud(obj_last);
% pc_new_last.Color=pc_last.Color;
% figure,pcshow(pc_new_last); xlabel('X'); ylabel('Y');zlabel('Z');



%% bottom
% roi = [  1 512  1   512  1 205 ]; %% 165 
% pc = select(pc,findPointsInROI(pc,roi));
% pcshow(pc)

roi = [  1 512  1   512  1 99 ]; %% 165 
pc_bottom = select(pc,findPointsInROI(pc,roi));
figure;pcshow(pc_bottom)


obj_bottom = im2double(pc_bottom.Location);
obj_bottom(:,3) =  obj_bottom(:,3)- 50;

pc_new_bottom = pointCloud(obj_bottom);
pc_new_bottom.Color=pc_bottom.Color;
figure,pcshow(pc_new_bottom); xlabel('X'); ylabel('Y');zlabel('Z');

%% merge all

ptCloud = pcmerge(pc_new,pcmerge(pc_new_bottom,pc_middle,1),1);
pcshow(ptCloud);

obj = im2double(ptCloud.Location);
obj_c = im2double(ptCloud.Color);


obj(:,1) = obj(:,1)/max(obj(:,1)) -0.5;
obj(:,2) = obj(:,2)/max(obj(:,2)) -0.5 ;
obj(:,3) = obj(:,3)/max(obj(:,3));

obj(:,1) = obj(:,1) /100;
obj(:,2) = obj(:,2) /100;

