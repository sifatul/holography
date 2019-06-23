%Main program of hologram generation by wavefront recording plane method
%by Openholo library project
%2017-10-30 update
%

addpath('D:\Mathlab\wrp\data') 
addpath('D:\Mathlab\wrp\data_scaled');
cd D:\Mathlab\wrp

clc;clear;close all;

 
% zhao
% obj(:,1) = (obj(:,1)/450);
% obj(:,2) = (obj(:,2)/450);  
% obj(:,3) = (obj(:,3)/20); 

 
three_object  
obj(:,1) = obj(:,1)*3.5; 
obj(:,2) = obj(:,2)*3.5;  
obj(:,3) = obj(:,3)*50;   

%% OBJECT = BEAR 

 %{
load ./data/bear_67w.txt;
obj_name = "bear_67w";

obj = zeros(size(bear_67w,1),3);
obj(:,1) = bear_67w(:,1);
obj(:,2) = bear_67w(:,2);
obj(:,3) = bear_67w(:,3);


sitha=-90;
betha=0;
gamma=90;

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


%figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

obj(:,1) = obj(:,1) /3000;
obj(:,2) = obj(:,2) /3000;
obj(:,3) = obj(:,3)/200;
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

% obj = ptCloud.Location;
% obj_c = im2double ( ptCloud.Color);
% obj = obj./3000 ;
% obj(:,2) = obj(:,2)./5;

% Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
% obj_X = obj(:,1);
% obj_y = obj(:,2);
% obj_z = obj(:,3);
% 
% 
% obj2(:,1) = obj_X(Visible_point_indexes);
% obj2(:,2) = obj_y(Visible_point_indexes);
% obj2(:,3) = obj_z(Visible_point_indexes);

ptCloud = pointCloud(obj); 
figure,pcshow(ptCloud); 
%}
 


%% OBJECT = Panda 
 
%{
load C_panda.txt;
 
obj_name = "C_panda";
% should try for d = 0.35; 

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

sitha=70;
betha=0;
gamma=90;

XRS=[1 0 0;
     0 cosd(sitha) -sind(sitha);
     0 sind(sitha) cosd(sitha);];
YRS=[cosd(betha) 0 -sind(betha) ;
    0 1 0 ;
    sind(betha) 0 cosd(betha) ;];
ZRS=[cosd(gamma) -sind(gamma) 0 ;
    sind(gamma) cosd(gamma) 0 ;
    0 0 1 ;]; 
obj2=obj*XRS;

  obj(:,1) = obj(:,1)./5 ;
  obj(:,2) = obj(:,2)./5 ; %less than 10
obj(:,3) = obj(:,3) ;


% Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
% obj_X = obj(:,1);
% obj_y = obj(:,2);
% obj_z = obj(:,3);
% obj2(:,1) = obj_X(Visible_point_indexes);
% obj2(:,2) = obj_y(Visible_point_indexes);
% obj2(:,3) = obj_z(Visible_point_indexes);

ptCloud2 = pointCloud(obj2); 
ptCloud2.Color= ptCloud.Color ;
figure,pcshow(ptCloud2);  xlabel('X');
ylabel('Y');
zlabel('Z');
%}
 

%% 
%{
 pc = pcread ('./data/box_ninja.ply');
obj = im2double(pc.Location);
obj_c = im2double(pc.Color);

obj(:,1) = (obj(:,1)+0.050)/100;
obj(:,2)= obj(:,2)/100;
obj(:,3)= obj(:,3)/100 ;

figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');
%}
%% 
%{ 
obj_name = "small_bear";
ptCloud = pcread ('bear3.ply');
obj = im2double( ptCloud.Location);
obj(:,1)=obj(:,1)+abs(max(obj(:,1))+min(obj(:,1)))/2-0.05;%+abs(min(obj1(:,1)))+1;
obj(:,2)=obj(:,2)+abs(max(obj(:,2))+min(obj(:,2)))/2-0.025;%+abs(min(obj1(:,2)))+1;
% obj(:,3)=obj(:,3)+abs(min(obj(:,3)))/2+1;
obj(:,1) = obj(:,1)/150;
obj(:,2) = obj(:,2)/150;
obj(:,3) = obj(:,3)/40;
ptCloud2 = pointCloud(obj);
ptCloud2.Color=uint8(ptCloud.Color);
figure,pcshow(ptCloud2);
xlabel('X');
ylabel('Y');
zlabel('Z');
obj_c = double(ptCloud2.Color)/256 ;
%}
 
%% 

 %{
d = 0.18; 
obj_name = "princes";
ptCloud = pcread ('./data/princes.ply');
obj = im2double( ptCloud.Location);
obj(:,1)=obj(:,1)+abs(max(obj(:,1))+min(obj(:,1)))/2-0.05;%+abs(min(obj1(:,1)))+1;
obj(:,2)=obj(:,2)+abs(max(obj(:,2))+min(obj(:,2)))/2 -0.1;%+abs(min(obj1(:,2)))+1;
% obj(:,3)=obj(:,3)+abs(min(obj(:,3)))/2+1;
obj(:,1) = obj(:,1)/80;
obj(:,2) = obj(:,2)/80;
obj(:,3) = obj(:,3)/10;
ptCloud2 = pointCloud(obj);
ptCloud2.Color=uint8(ptCloud.Color);
figure,pcshow(ptCloud2);
xlabel('X');
ylabel('Y');
zlabel('Z');
obj_c = double(ptCloud2.Color)/256 ;
%}
%% 


obj_depth = max(obj) - min(obj)

 
file_type = '.bmp';
project_name = mfilename;

Hologram_resolution_x = 1980;  
Hologram_resolution_y = 1024;  % Hologram resolution  
Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;
%Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
Hologram_sampling_interval = 7.4e-6; %8.4e-6;%7.4e-6;8.4e-6; %;            % Hologram sampling interval




% Hologram_wrp  =gpuArray( zeros(Hologram_resolution, Hologram_resolution, 3));  %gpuArray( zeros(Hologram_resolution));

lambda = [632.8e-9 532e-9 473e-9]; %RGB; 
 

obj_z = obj(:,3);
%  obj_z_temp = obj(:,3);
% obj_z= round(obj_z_temp,3);
Cut = sort(unique(obj_z));

ROWS= Hologram_resolution_x;                                     
COLS= Hologram_resolution_y;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
obj_depth= max(obj(:,3)) -min(obj(:,3));
d = 0.5; 

Nx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);  
Ny = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );
counter_limit = [50]; 
max_counter= 0;


sub_dir = strcat(project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)));
% % mkdir(sub_dir );
% focus_location ='face_';
% file_name = strcat(focus_location,'time_',sprintf('%d',Hologram_resolution),'.txt');
% full_file_name = fullfile(sub_dir, file_name);
% fileID = fopen(full_file_name,'w'); 

for ii=1:length(counter_limit)
   
counter2_limit = counter_limit(ii);
original=zeros(Hologram_resolution_x,Hologram_resolution_y,3); 
original2=zeros(Hologram_resolution_x,Hologram_resolution_y,3); 
 tic
for color_index = 1:length(lambda)
%     color_index
    counter = 0;
    k = 2*pi/lambda(color_index); 
    WRP =  zeros(Hologram_resolution_x,Hologram_resolution_y);
    WRPHologram =   zeros(Hologram_resolution_x,Hologram_resolution_y);
    z_rem = -1; 
    total_wrp =0;
%     [a,b]=hist(obj_z,unique(obj_z));
%     xlabel('Depth Layer');
% ylabel('');


    counter2 =0; 
  
    for o = 1: length(Cut)
        indexes = find(obj_z==Cut(o));
        if length(obj_z(indexes)) ~= z_rem  
            z_rem = length(obj_z(indexes));              
%             total_wrp = total_wrp +1;    
            z_wrp = Cut(o)-0.0005;  %  WRP location   
             if(o~=1)
                d_wrp_to_hologram = max(abs(d-z_wrp));
%                 d_wrp_to_hologram
                WRPHologram = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
                WRP =  zeros(Hologram_resolution_x,Hologram_resolution_y);
                  
             end           
                        
%             if(counter2>max_counter)
%                 max_counter  = counter2;
%             end
%             counter2 = 0;             
        end 
%         counter2 = counter2 +1;
        z =  (z_wrp -  Cut(o));
        color_for_this_depth = obj_c(indexes,color_index);

        active_area =  (round(abs(lambda(color_index).*z(1)./(Hologram_sampling_interval^2)/2)+0.5).*2-1 );        %sampling size of N 
        Nx_for_this_depth = Nx(indexes);
        Ny_for_this_depth = Ny(indexes);

        [y_run, x_run]= meshgrid((-(active_area-1)/2:(active_area-1)/2)*Hologram_sampling_interval,(-(active_area-1)/2:(active_area-1)/2)*Hologram_sampling_interval);
        r =  (sign(z(1))*sqrt(z(1)^2 + y_run.^2 + x_run.^2)); 
        for Ni = 1 : length( indexes)
            nx = Nx_for_this_depth(Ni); ny = Ny_for_this_depth(Ni);  
            counter = counter +1;
            fprintf('counter: %d color: %d\n',counter ,color_index);           
            Sub_hologram =  color_for_this_depth(Ni)*exp(1j*rand*2*pi)*exp(1j*k*r)./r;%* exp(1j*rand*2*pi)
            WRP(nx:nx+active_area-1,ny:ny+active_area-1)= WRP(nx:nx+active_area-1,ny:ny+active_area-1) + Sub_hologram;
         
        end 
%            max(max(WRP))
   
    end %end of all layers
    
 WRPHologram = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
    
% end_time= toc;        
   close all    
 for p = 31%reconstructed    
    d2 = d+0.002 - p*0.001; 
    original_red =  FresnelPropogation(k,v, h,-d2,WRPHologram);
    original_abs_red= abs((original_red)); 
    original(:,:,color_index)=original_abs_red;
    original2(:,:,color_index)= 255.*(original_abs_red./max(max(original_abs_red))) ;
%      figure; imshow(original_abs_red,[]);  
%      
%      file_name = strcat(focus_location,project_name,'_color_',num2str(color_index),'_d_',num2str(d),'total_wrp_',num2str(total_wrp),'_recon_d_',num2str(d2),file_type);
%     fullFileName = fullfile(sub_dir,file_name); 
%     save_hologram(gather(WRPHologram),fullFileName);
% %      
end
    
%   disp( color_index);

%     title(num2str(color_index));
end
 
%% reconstruction for red

end_time= toc;
imshow(rot90(uint8((original2 ))),[]);


% fprintf(fileID,' %d %0.1f \n ',total_wrp, end_time);
% file_name = strcat('full_hologram',project_name,'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
% fullFileName = fullfile(sub_dir,file_name); 
% save_hologram(gather(WRPHologram),fullFileName);
% % imshow(rot90(uint8(gather(original))),[]);
file_name = strcat(focus_location,project_name,'_d_',num2str(d),'total_wrp_',num2str(total_wrp),'_recon_d_',num2str(d2),file_type);
fullFileName = fullfile(sub_dir,file_name); 
imwrite((uint8((original2))),fullFileName);

fprintf(fileID,' %d %0.1f \n ',total_wrp, end_time);
%  disp(counter2_limit);
end
 
fclose(fileID);
