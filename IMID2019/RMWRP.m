%Main program of hologram generation by wavefront recording plane method
%by Openholo library project
%2017-10-30 update
%


clc;clear;close all;
 %{
pcread ('bunny.ply');
obj_name = "bunny";
  temp = ans.Location;
    obj(:,1) = temp(:,2);
    obj(:,2) = temp(:,1);
    obj(:,3) = temp(:,3);
    obj(:,2) = obj(:,2) +0.05;
    obj(:,1) = obj(:,1) -0.1;
    obj = (obj )/80;
 figure; plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

 %}

 %{
  pcread ('horse.ply');
  obj_name = "horse";
  temp = ans.Location;
    obj(:,1) = temp(:,2);
    obj(:,2) = temp(:,3);
    obj(:,3) = temp(:,1);
 figure; plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');
 obj1 =  unique(obj,'rows');
  figure; plot3(obj1(:,1),obj1(:,2),obj1(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');
 obj = (obj+0.05)./80;
 
 %}

%{
load tea_pot.mat;
obj_name = "tea_pot";


obj(:,1)=(Obj(:,1)/40000 -0.0002 ).*1;
obj(:,2)= (Obj(:,2)/20000 )+1.5;
obj(:,3)= (Obj(:,3)/30000);

figure; plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

sitha=180;
betha=0;
gamma=180;

XRS=[1 0 0;
     0 cosd(sitha) -sind(sitha);
     0 sind(sitha) cosd(sitha);];
YRS=[cosd(betha) 0 -sind(betha) ;
    0 1 0 ;
    sind(betha) 0 cosd(betha) ;];
ZRS=[cosd(gamma) -sind(gamma) 0 ;
    sind(gamma) cosd(gamma) 0 ;
    0 0 1 ;]; 
obj=obj*XRS;

obj(:,1)=obj(:,1)+abs(max(obj(:,1))+min(obj(:,1)))/2+1;%+abs(min(obj1(:,1)))+1;
obj(:,2)=obj(:,2)+abs(max(obj(:,2))+min(obj(:,2)))/2+1;%+abs(min(obj1(:,2)))+1;
obj(:,3)=obj(:,3)+abs(min(obj(:,3)))/2+1;

ptCloud=pointCloud(obj);
% ptCloud.Color=uint8(obj_c);
figure,pcshow(ptCloud);
xlabel('X');
ylabel('Y');
zlabel('Z');
%}
%% 
% min_z = min(obj(:,2));
% max_z = max(obj(:,2));
% mid_z = (max_z-min_z)./2;
% 
% obj1= obj(obj(:,2)<= min_z+mid_z,:); 

%figure; plot3(obj1(:,1),obj1(:,2),obj1(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');
 
%load obj;
%{
load dragon_0;
obj_name = "dragon_0";
obj = (ptCloud2.Location)./10 ;
obj(:,1) = (obj(:,1))./6 ;
obj(:,2) = (obj(:,2))./15 ;
obj(:,3) = (obj(:,3)).*8 ;
%}
%{
tic

img1=rgb2gray(imread('lena256_color.bmp'));
img=double(imread('lena256_color.bmp'));

obj=depthimg2point(rgb2gray(img),0);

obj_x=obj(:,1);
obj_y=obj(:,2);
obj_z=obj(:,3);
% obj(length(obj_z)/2:end,3)=obj(length(obj_z)/2:end,3).*2;
obj_c=zeros(length(obj_x),3);

for a=1:length(obj(:,1))
    obj_c(a,:)=reshape(img(obj_x(a),obj_y(a),:),[],3);
end

figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);

figure,pcshow(ptCloud);

obj_min=min(obj);
obj_center=(max(obj)+min(obj))./2;
obj(:,1)= (obj(:,1)/max(obj(:,1)))-0.4;
obj(:,2)= obj(:,2)/max(obj(:,2))-0.4;
%obj(:,3)= obj(:,3)/max(obj(:,3)); 
obj(:,3) =  rand(65536,1) ;
%obj= obj.*2.4; 
obj = obj./1000 ;
obj(:,3) = obj(:,3).*1000 ;
obj(:,3) = round(obj(:,3).*100);
%lambda = 532e-9; %green wave  
%}
%% OBJECT = BEAR 
 %{

load bear_67w.txt;
obj_name = "bear_67w";

obj = zeros(size(bear_67w,1),3);
obj(:,1) = bear_67w(:,1);
obj(:,2) = bear_67w(:,2);
obj(:,3) = bear_67w(:,3);
%figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


obj_c=zeros(length(obj),3);

obj_c(:,1) = bear_67w(:,4).*256;
obj_c(:,2) = bear_67w(:,5).*256;
obj_c(:,3) = bear_67w(:,6).*256; 
ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);
figure,pcshow(ptCloud); 

xlabel('X');
ylabel('Y');
zlabel('Z');

obj = ptCloud.Location;
obj_c = im2double ( ptCloud.Color);
obj = obj./3000 ;
obj(:,2) = obj(:,2)./5;

Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
obj_X = obj(:,1);
obj_y = obj(:,2);
obj_z = obj(:,3);
obj2(:,1) = obj_X(Visible_point_indexes);
obj2(:,2) = obj_y(Visible_point_indexes);
obj2(:,3) = obj_z(Visible_point_indexes);

ptCloud = pointCloud(obj2); 
figure,pcshow(ptCloud); 
%}
 


%% OBJECT = Panda 
 %{
load C_panda.txt;
 
obj_name = "C_panda";


obj = zeros(size(C_panda,1),3);
obj(:,1) = C_panda(:,1);
obj(:,2) = C_panda(:,2);
obj(:,3) = C_panda(:,3);
%figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


obj_c=zeros(length(obj),3);

obj_c(:,1) = C_panda(:,5);
obj_c(:,2) = C_panda(:,6);
obj_c(:,3) = C_panda(:,7);
ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);
figure,pcshow(ptCloud); 

xlabel('X');
ylabel('Y');
zlabel('Z');

obj = ptCloud.Location;
obj_c = im2double ( ptCloud.Color);
 obj = obj./1000 ;

Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
obj_X = obj(:,1);
obj_y = obj(:,2);
obj_z = obj(:,3);
obj2(:,1) = obj_X(Visible_point_indexes);
obj2(:,2) = obj_y(Visible_point_indexes);
obj2(:,3) = obj_z(Visible_point_indexes);

ptCloud = pointCloud(obj2); 
figure,pcshow(ptCloud); 
%}
%% %% OBJECT = three_obj


  

load three_obj.txt

obj_name = "three_obj";


obj = zeros(size(three_obj,1),3);
obj(:,1) = three_obj(:,1);
obj(:,2) = three_obj(:,2);
obj(:,3) = three_obj(:,3);

obj(:,1) = (obj(:,1)/max(obj(:,1)))+0.25;
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3));

obj = obj/1000; 
obj = obj*1.30;
%figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


obj_c=zeros(length(obj),3);

obj_c(:,1) = three_obj(:,4).*256;
obj_c(:,2) = three_obj(:,5).*256;
obj_c(:,3) = three_obj(:,6).*256; 
ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);
figure,pcshow(ptCloud); 



xlabel('X');
ylabel('Y');
zlabel('Z');

obj = ptCloud.Location;
obj_c = im2double ( ptCloud.Color);

%}

%% OBJECT = LENA
%{
 obj_name = "LENA";
% img1=rgb2gray(imread('lena256_color.bmp'));
img=double(imread('lena256_color.bmp'));

obj=depthimg2point(rgb2gray(img),0);

obj_x=obj(:,1);
obj_y=obj(:,2);
obj_z=obj(:,3);
% obj(length(obj_z)/2:end,3)=obj(length(obj_z)/2:end,3).*2;
obj_c=zeros(length(obj_x),3);

for a=1:length(obj(:,1))
    obj_c(a,:)=reshape(img(obj_x(a),obj_y(a),:),[],3);
end

figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);

figure,pcshow(ptCloud);

obj_min=min(obj);
obj_center=(max(obj)+min(obj))./2;
obj(:,1)= (obj(:,1)/max(obj(:,1)))-0.4;
obj(:,2)= obj(:,2)/max(obj(:,2))-0.4;
%  obj(:,3)= obj(:,3)/max(obj(:,3)); 
obj(:,3) =  rand(65536,1) ;
% %obj= obj.*2.4; 
% obj = obj./1000 ;

%obj(:,3) =  rand(65536,1) ;
%obj= obj.*2.4; 
obj = (obj-0.2)./400 ;
 
%}

%% 



tic
file_type = '.bmp';
project_name='RMWRP';
tic
Hologram_resolution = 512;                       % Hologram resolution    
% Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
Hologram_sampling_interval = 7.4e-6;            % Hologram sampling interval



% Hologram_wrp  =gpuArray( zeros(Hologram_resolution, Hologram_resolution, 3));  %gpuArray( zeros(Hologram_resolution));

lambda = [632.8e-9 532e-9 473e-9]; %RGB; 
 

obj_z = obj(:,3);
 
Cut = sort(unique(obj_z));
original=zeros(Hologram_resolution,Hologram_resolution,3); 

ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));


d=0.06;
d2 = d+0.002 - 7*0.0002;

Nx = gpuArray( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution)/2);  
Ny = gpuArray( round(obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution)/2 );

for color_index = 1:length(lambda)
    k = 2*pi/lambda(color_index); 
    Hologram_wrp = gpuArray( zeros(Hologram_resolution));
    z_rem = -1; 
    counter =0; 
    for o = 1: length(Cut)
        indexes = find(obj_z==Cut(o));
        if length(obj_z(indexes)) ~= z_rem
            z_rem = length(obj_z(indexes));
            z_wrp = obj_z(indexes)+0.0005;  %  WRP location        
            counter = counter +1;
        end 
        z = gpuArray (z_wrp -  Cut(o));
        color_for_this_depth = obj_c(indexes,color_index);

        N = gpuArray (round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 );        %sampling size of N 
        Nx_for_this_depth = Nx(indexes);
        Ny_for_this_depth = Ny(indexes);

        for Ni = 1 : length(N)
            nx = Nx_for_this_depth(Ni); ny = Ny_for_this_depth(Ni);
            [y_run, x_run]= meshgrid((-(N(Ni)-1)/2:(N(Ni)-1)/2)*Hologram_sampling_interval,(-(N(Ni)-1)/2:(N(Ni)-1)/2)*Hologram_sampling_interval);
            r = gpuArray (sign(z(Ni))*sqrt(z(Ni)^2 + y_run.^2 + x_run.^2));
            Sub_hologram =  color_for_this_depth(Ni)* exp(1j*2*pi)*exp(1j*k*r)./r;
            Hologram_wrp(nx:nx+N(Ni)-1,ny:ny+N(Ni)-1) = Hologram_wrp(nx:nx+N(Ni)-1,ny:ny+N(Ni)-1)+ Sub_hologram;
        end

    end


    
    
    WRPHologram  = FresnelPropagation(gather(Hologram_wrp), Hologram_sampling_interval, Hologram_sampling_interval, d, lambda(color_index));
    original_red =  FresnelPropogation(k,v, h,-d2,WRPHologram);
    original_abs_red= abs(original_red); 
    original(:,:,color_index)=255.*(original_abs_red./max(max(original_abs_red)));
    figure; imshow(rot90(abs(original(:,:,color_index))),[]);
 
    title(['image o = ',num2str(color_index)])
     
 
    
end
 
%% reconstruction for red
  


toc
% file_prefix = 'hologram_red';
% file_name = strcat(obj_name,project_name,sprintf('%0.1f',Hologram_resolution),sprintf('%0.9f',Hologram_sampling_interval),file_prefix,file_type);
% save_hologram(WRPHologram_red,file_name);
imshow(rot90(uint8(gather(original))),[]);
 file_prefix = 'full_color_reconstruction';
file_name = strcat(obj_name,project_name,sprintf('%0.1f',Hologram_resolution),sprintf('%0.9f',Hologram_sampling_interval),file_prefix,file_type);
imwrite(rot90(uint8(gather(original))),file_name);
