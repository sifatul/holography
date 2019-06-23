addpath('D:\Mathlab\wrp\data'); 
mkdir('D:\Mathlab\wrp\data')
depth=double(imread('papillon_depth.png'));
color=double(imread('papillon_color.jpg'));
obj=depthimg2point(depth,0);

obj_x=obj(:,1);
obj_y=obj(:,2); 
% obj(1:end,3)= rand(length(obj_x),1 ) ;
obj_c=zeros(length(obj_x),3);

for a=1:length(obj(:,1))
    obj_c(a,:)=reshape(color(obj_x(a),obj_y(a),:),[],3);
end

ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);

 

figure,pcshow(ptCloud);
xlabel('X');
ylabel('Y');
zlabel('Z');
pcwrite(ptCloud,'papilon2.ply');

% obj_z= obj(:,3);
% cut = unique(obj_z);
%  
% [a,b]=hist(cut,unique(cut))
% 
% obj_min=min(obj);
% obj_center=(max(obj)+min(obj))./2;
% obj(:,1)= (obj(:,1)/max(obj(:,1)))-0.4;
% obj(:,2)= obj(:,2)/max(obj(:,2))-0.4;
% %  obj(:,3)= obj(:,3)/max(obj(:,3)); 
% % obj(:,3) = 
% % %obj= obj.*2.4; 
% % obj = obj./1000 ;
% 
% %obj(:,3) =  rand(65536,1) ;
% %obj= obj.*2.4; 
% obj = (obj-0.2)./400 ;
 