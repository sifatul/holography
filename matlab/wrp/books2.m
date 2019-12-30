 
clc;clear all;close all;
data = load ('C:\Users\user\Desktop\new_data\books2.txt');
obj_name = strcat('books2');

% data = data_load(data_load(:,3)>130 &  data_load(:,2)>-50,:);
obj2(:,1) = data(:,1);
obj2(:,2) = data(:,2);
obj2(:,3) = data(:,3);

obj2_c=zeros(length(obj2),3);
obj2_c(:,1) = data(:,4).*256;
obj2_c(:,2) = data(:,5).*256;
obj2_c(:,3) = data(:,6).*256; 

del_idx = find(obj2(:,3)>295& obj2(:,3)<370);
obj2(del_idx,:)= [];
obj2_c(del_idx,:) = [];

temp_cut = unique(obj2(:,3));
idx = find( obj2(:,3) == max(temp_cut));
obj2(idx,:) = [];
obj2_c(idx,:)=[];
temp_cut = unique(obj2(:,3));
idx = find( obj2(:,3) == max(temp_cut));
obj2(idx,:) = [];
obj2_c(idx,:)=[];



pc = pointCloud(obj2);
pc.Color=uint8(obj2_c);
figure,pcshow(pc); xlabel('X'); ylabel('Y');zlabel('Z');

roi = [    NaN   98  -80 30  NaN 417 ]; %% 165 

indices = findPointsInROI(pc,roi); 

pc = select(pc,indices);
figure,pcshow(pc); xlabel('X'); ylabel('Y');zlabel('Z');


obj = im2double(pc.Location);

obj(:,1) = (obj(:,1)/max(obj(:,1)));
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3)); 

obj_c= im2double(pc.Color).*256;
unique(obj(:,3))
