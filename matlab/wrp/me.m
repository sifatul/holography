 close all; clear all; clc;

temp = load('C:\Users\user\Desktop\new_data\sahinur.txt '); 
obj_name = "me";
% 
% obj_temp = zeros(size(temp,1),3);
% obj_temp(:,1) = temp(:,1);
% obj_temp(:,2) = temp(:,2);
% obj_temp(:,3) = temp(:,3);
% 
% obj_temp_c=zeros(length(obj_temp),3);
% 
% obj_temp_c(:,1) = temp(:,4).*256;
% obj_temp_c(:,2) = temp(:,5).*256;
% obj_temp_c(:,3) = temp(:,6).*256;  
% ptCloud = pointCloud(obj_temp );
% ptCloud.Color=uint8(obj_temp_c );
% figure,pcshow(ptCloud); 
% 
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% 
% roi = [  -400 600  -300   100  400 1200 ]; %% 165
%  
% indices = findPointsInROI(ptCloud,roi); 
% 
% pc = select(ptCloud,indices);
% % 
% % 
% % 
% figure,pcshow(pc); xlabel('X'); ylabel('Y');zlabel('Z');
% obj = im2double(pc.Location);  
% obj_c = im2double(pc.Color).*256;
% idx_rem = find(obj(:,3)>=800 &  obj(:,1) < 50 |  obj(:,3)<=500);
% obj(idx_rem,:) = [];
% obj_c(idx_rem,:) = [];
% % 
% % 
% pc2 = pointCloud(obj );
% pc2.Color=uint8(obj_c );
% figure,pcshow(pc2);  
% 
% 
% 
%  idx = im2double(pc2.Location);  
% obj_c = im2double(pc2.Color).*256; 
% 
% 
% idx = find(obj(:,3)>500 & obj(:,3)<580 );
% obj(idx,3) = max(obj(idx,3) );
% % obj(idx,3) = round((obj(:,3)/max(obj(:,3))),4);
% 
% 
% % obj(:,1) = (obj(:,1)/max(obj(:,1))); 
% % obj(:,2) = obj(:,2)/max(obj(:,2));
% % obj(:,3) = obj(:,3)/max(obj(:,3));
% 
% 
% % obj(:,1) = (obj(:,1)/100);
% % obj(:,2) = (obj(:,2)/100); 
% % obj(:,3) = (obj(:,3) );
% 
% pc2 = pointCloud(obj );
% pc2.Color=uint8(obj_c );
% figure,pcshow(pc2);  
pc= pcread('me_sahinur_increase.ply');
% figure; pcshow(pc);
obj= im2double(pc.Location);
obj_c= im2double(pc.Color);
idx_rem = find(obj(:,3)>=800 &  obj(:,1) < 50 |  obj(:,2) < -200 );
obj(idx_rem,:) = [];
obj_c(idx_rem,:) = [];

%delete 

 idx_rem = find(obj(:,3)<580) ;
 obj(idx_rem,3) = max(obj(idx_rem,3));
 
  idx_rem = find(obj(:,3)>800 & obj(:,3)<960) ;
 obj(idx_rem,3) = max(obj(idx_rem,3));
 
 
 idx_rem = find(obj(:,3)>=max(obj(:,3)) );
obj(idx_rem,:) = [];
obj_c(idx_rem,:) = [];

 
pc2 = pointCloud(obj );
pc2.Color=uint8(obj_c*256 );
figure,pcshow(pc2);  
xlabel('X');
ylabel('Y');
zlabel('Z');
% 

obj(:,1) = (obj(:,1)/max(obj(:,1))); 
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3));


obj(:,1) = (obj(:,1)/100);
obj(:,2) = (obj(:,2)/100); 
obj(:,3) = (round(obj(:,3),3) );
 idx_rem = find(obj(:,3)>=max(obj(:,3)) );
obj(idx_rem,:) = [];
obj_c(idx_rem,:) = [];

unique(obj(:,3))
pc2 = pointCloud(obj );
pc2.Color=uint8(obj_c*256 );
figure,pcshow(pc2);  
%  
