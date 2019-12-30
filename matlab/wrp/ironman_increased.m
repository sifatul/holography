 
clc;clear all;close all; 
addpath('D:\Mathlab\wrp\data_scaled');
obj_name = strcat('ironman_increased');
 pc = pcread('ironman_increased.ply');

obj = im2double(pc.Location);

% half_idx = find(obj(:,3)<460);
% obj(half_idx,3) = obj(half_idx,3) - 150;

obj(:,1) = (obj(:,1)/max(obj(:,1)));
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3)); 

obj_c = im2double(pc.Color).*256;
% idx = find( obj(:,3) == max(obj(:,3)));
% obj(idx,:) = [];
% obj_c(idx,:)=[];

obj(:,3) = round( obj(:,3),3);

unique(obj(:,3))

ptCloud2 = pointCloud(obj);
ptCloud2.Color=uint8(obj_c);
figure,pcshow(ptCloud2); xlabel('X'); ylabel('Y');zlabel('Z');
