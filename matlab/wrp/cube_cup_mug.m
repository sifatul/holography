 clc; clear all; close all;
data = load ('C:\Users\user\Desktop\new_data\cup_mug_cube3.txt');
obj_name = strcat('cup_cube_mug');

% data = data_load(data_load(:,3)>130 &  data_load(:,2)>-50,:);
obj2(:,1) = data(:,1);
obj2(:,2) = data(:,2);
obj2(:,3) = data(:,3);

obj2_c=zeros(length(obj2),3);
obj2_c(:,1) = data(:,4).*256;
obj2_c(:,2) = data(:,5).*256;
obj2_c(:,3) = data(:,6).*256; 

% temp_cut = unique(obj2(:,3));
% idx = find( obj2(:,3) == min(temp_cut));
% obj2(idx,:) = [];
% obj2_c(idx,:)=[];
% temp_cut = unique(obj2(:,3));
% idx = find( obj2(:,3) == max(temp_cut));
% obj2(idx,:) = [];
% obj2_c(idx,:)=[];



ptCloud = pointCloud(obj2);
ptCloud.Color=uint8(obj2_c);
figure,pcshow(ptCloud); xlabel('X'); ylabel('Y');zlabel('Z');

 roi = [    -145.1828 98 -30 68  151 280 ]; %% 165



indices = findPointsInROI(ptCloud,roi); 

pc = select(ptCloud,indices);

figure,pcshow(pc); xlabel('X'); ylabel('Y');zlabel('Z');


obj = im2double(pc.Location);

obj(:,1) = (obj(:,1)/max(obj(:,1)));
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3)); 

obj_c= im2double(pc.Color).*256;
 
[Cut,~,idx] = unique(obj(:,3));
