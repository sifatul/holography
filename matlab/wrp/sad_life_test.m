addpath('D:\Mathlab\wrp\data')  
close all;
% 
% obj_name = 'try1_100';
%   
% pc = pcread('try1_100.ply');
% 
% 
% obj = im2double(pc.Location);
% 
% obj(:,1) = (obj(:,1)/max(obj(:,1)));
% obj(:,2) = obj(:,2)/max(obj(:,2));
% obj(:,3) = obj(:,3)/max(obj(:,3)); 
% 
% obj_c= im2double(pc.Color).*256;
% 
% figure,pcshow(pc);  xlabel('X'); ylabel('Y');zlabel('Z');
% 
% 
% [temp_cut,~,idx] = unique(obj(:,3));
% n = accumarray(idx(:),1);
% index = find(n<44);
% range = temp_cut(index);
% Lia = ismember(obj(:,3),range);
% idx = find(Lia);
% obj(idx,:) = [];
% obj_c(idx,:)=[];
% 
% 
% [temp_cut,~,idx] = unique(obj(:,3));
% n = accumarray(idx(:),1);
% index = find(n==46);
% range = temp_cut(index(1));
% Lia = ismember(obj(:,3),range);
% idx = find(Lia);
% obj(idx,:) = [];
% obj_c(idx,:)=[]; 
% 
% 
% ptCloud2 = pointCloud(obj);
% ptCloud2.Color=uint8(obj_c);
% figure,pcshow(ptCloud2); xlabel('X'); ylabel('Y');zlabel('Z');
% 


% addpath('D:\Mathlab\wrp\data_scaled');
% ptCloudA = pcread('book.ply');
% ptCloudB = pcread('cube.ply');
% ptCloudC = pcread('cup.ply');
% % ptCloudD = pcread('chicken.ply');
% 
% ptCloud = pcmerge(ptCloudA,ptCloudC,1);
% figure,pcshow(ptCloud); 
% xlabel('X'); ylabel('Y');zlabel('Z');
% 
% 
% obj_cube = im2double(ptCloudC.Location);
% obj_cube(:,2) = obj_cube(:,2)-40;
% 
% ptCloud_cube = pointCloud(obj_cube);
% ptCloud_cube.Color=(ptCloudC.Color);
% figure,pcshow(ptCloud_cube); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% ptCloud = pcmerge(ptCloudA,ptCloudB,1);
% figure,pcshow(ptCloud); 
% 
% 
% cube = im2double(ptCloudB.Location);
% 
% cube(:,1) = cube(:,1)+30 ;
% cube(:,2) = cube(:,2) ;
% cube(:,3) = cube(:,3)+350 ;
% 
% ptCloud_cube = pointCloud(cube);
% ptCloud_cube.Color=(ptCloudB.Color);
% 
% ptCloud2 = pcmerge(ptCloudA,ptCloud_cube,1);
% figure,pcshow(ptCloud2); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% cup_loc = im2double(ptCloudC.Location);
% cup_loc(:,3) = cup_loc(:,3) + 250;
% cup_loc(:,2) = cup_loc(:,2) -30;
% pt_cup = pointCloud( cup_loc);
% 
% pt_cup.Color=(ptCloudC.Color);
% 
% ptCloud2 = pcmerge(ptCloud2,pt_cup,1);
% figure,pcshow(ptCloud2); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% 
% 
% 
% 
% 
% 
% 
% chick_loc = im2double(ptCloudD.Location);
% chick_loc(:,1) = chick_loc(:,1)-80;
% chick_loc(:,2) = chick_loc(:,2)-5;
% chick_loc(:,3) = chick_loc(:,3)+50;
% ptCloud_chick = pointCloud(chick_loc);
% ptCloud_chick.Color=(ptCloudD.Color);
% 
% ptCloud3 = pcmerge(ptCloud2,ptCloud_chick,1);
% figure,pcshow(ptCloud3); xlabel('X'); ylabel('Y');zlabel('Z');
% 
% xlabel('X'); ylabel('Y');zlabel('Z');
% obj = im2double(ptCloud3.Location);
% Cut = unique(obj(:,3));
% figure,scatter([1:length(Cut)],Cut);
% address = 'D:\Mathlab\wrp\data_scaled\try1_0';
% pcwrite(ptCloud2,address);

close all;


addpath('D:\Mathlab\wrp\data_scaled');
ptCloudA = pcread('book.ply');
ptCloudB = pcread('cube.ply');
ptCloudC = pcread('cup.ply');
addpath('D:\Mathlab\libs') ;
% ptCloudD = pcread('chicken.ply');


  
obj_name = 'sad_life_test';



book = im2double(ptCloudA.Location);

 
book(:,3) =0;
ptCloud_book = pointCloud(book);
ptCloud_book.Color=(ptCloudA.Color);


cube = im2double(ptCloudB.Location);

cube(:,1) = cube(:,1)+30;
cube(:,2) = cube(:,2) ;
cube(:,3) = cube(:,3);

ptCloud_cube = pointCloud(cube);
ptCloud_cube.Color=(ptCloudB.Color);

ptCloud2 = pcmerge(ptCloud_book,ptCloud_cube,1);
figure,pcshow(ptCloud2); xlabel('X'); ylabel('Y');zlabel('Z');

cup_loc = im2double(ptCloudC.Location);


cup_loc(:,3) = cup_loc(:,3);
cup_loc(:,2) = cup_loc(:,2) -30;
pt_cup = pointCloud( cup_loc);

pt_cup.Color=(ptCloudC.Color);

pc = pcmerge(ptCloud2,pt_cup,1);
figure,pcshow(pc); xlabel('X'); ylabel('Y');zlabel('Z');



  
% pc = pcread('try1_100.ply');

obj = im2double(pc.Location); 
obj_c= im2double(pc.Color).*256;



[temp_cut,~,idx] = unique(obj(:,3));
n = accumarray(idx(:),1);
index = find(n<44);
range = temp_cut(index);
Lia = ismember(obj(:,3),range);
idx = find(Lia);
obj(idx,:) = [];
obj_c(idx,:)=[];


[temp_cut,~,idx] = unique(obj(:,3));
n = accumarray(idx(:),1);
% index = find(n==46);
% range = temp_cut(index(1));
% Lia = ismember(obj(:,3),range);
% idx = find(Lia);
% obj(idx,:) = [];
% obj_c(idx,:)=[]; 

 pc = pointCloud( obj);

pc.Color=(uint8(obj_c));

figure,pcshow(pc);  xlabel('X'); ylabel('Y');zlabel('Z');
% pcwrite(pc,'D:\Mathlab\wrp\data_scaled\sad_life.ply');

% close all;
obj = im2double(pc.Location);

obj(:,1) = (obj(:,1)/max(obj(:,1)));
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3)); 

obj_c= im2double(pc.Color).*256;

length(unique(obj(:,3)))
% ptCloud2 = pointCloud(obj);
% ptCloud2.Color=uint8(obj_c);
% figure,pcshow(ptCloud2); xlabel('X'); ylabel('Y');zlabel('Z');

% 
% obj_parts = round((max(cut2)-min(cut2))/3);
% temp_cut = cut2;