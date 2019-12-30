close all; close all; clc;

addpath('D:\Mathlab\wrp\data_scaled'); 
addpath('D:\Mathlab\libs') ;
obj_name = 'temp';
   
% d = 0.23;
four_items = load ('four_items.txt');
obj2(:,1) = four_items(:,1);
obj2(:,2) = four_items(:,2);
obj2(:,3) = four_items(:,3);
 
obj2_c(:,1) = four_items(:,4).*256;
obj2_c(:,2) = four_items(:,5).*256;
obj2_c(:,3) = four_items(:,6).*256; 



[temp_cut,~,idx] = unique(obj2(:,3));
n = accumarray(idx(:),1);
index = find(n<30);
range = temp_cut(index);
Lia = ismember(obj2(:,3),range);
idx = find(Lia);
obj2(idx,:) = [];
obj2_c(idx,:)=[];
[temp_cut,~,idx] = unique(obj2(:,3));
temp_idx= find(obj2(:,3)==temp_cut(end));
obj2(temp_idx,:) = [];
obj2_c(temp_idx,:)=[];



ptCloud = pointCloud(obj2);
ptCloud.Color=uint8(obj2_c);
figure,pcshow(ptCloud); xlabel('X'); ylabel('Y');zlabel('Z');

% roi = [  -111.7362 167.3661  -44   46.2852   98 446.0000 ]; %% 123

% roi = [  -111.7362 167.3661  -60   46.2852  100 446.0000 ];  %%155 layers


% roi = [  -111.7362 167.3661  -50   46.2852  100 446.0000 ]; %%144


 roi = [  -111.7362 167.3661  -70   46.2852  102 446.0000 ]; %% 165



indices = findPointsInROI(ptCloud,roi); 

temp_pc = select(ptCloud,indices);



figure,pcshow(temp_pc); xlabel('X'); ylabel('Y');zlabel('Z');


% roi_cube = [  -111.7362 167.3661  -50   46.2852  350 446.0000 ]; %%144
% 
% indices_cube = findPointsInROI(ptCloud,roi_cube); 
%  
% 
% pc_cube = select(ptCloud,indices_cube);
% figure,pcshow(pc_cube); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% 
% roi_book = [  -111.7362 167.3661  -50   46.2852  200 350 ]; %%144
% indices_book = findPointsInROI(ptCloud,roi_book); 
% pc_book = select(ptCloud,indices_book);
% figure,pcshow(pc_book); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% 
% 
% roi_cup = [  -111.7362 167.3661  -50   46.2852  120 200 ]; %%144
% 
% indices_cup = findPointsInROI(ptCloud,roi_cup); 
% 
% pc_cup = select(ptCloud,indices_cup);
% 
% cup = im2double(pc_cup.Location);
% 
% cup(:,3) = cup(:,3); 
% pccup_new = pointCloud( cup);
% 
% pccup_new.Color=(pc_cup.Color);
%  
% 
% figure,pcshow(pccup_new); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% 
% 
% 
% book = im2double(pc_book.Location);
% 
% book(:,3) = book(:,3); 
% pc_book_new = pointCloud( book);
% 
% pc_book_new.Color=(pc_book.Color);
% 
% pc1 = pcmerge(pc_cube,pc_book_new,1);
% pc = pcmerge(pc1,pccup_new,1);
figure,pcshow(temp_pc); xlabel('X'); ylabel('Y');zlabel('Z');
obj = im2double(temp_pc.Location);


obj(:,1) = (obj(:,1)/max(obj(:,1)));
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3)); 

obj_c= im2double(temp_pc.Color).*256;

Cut= unique(obj(:,3));