addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
clc;clear;close all;

% papilon
% obj(:,1) = obj(:,1)/2.25;
% obj(:,2) = obj(:,2)/2.25;
% obj(:,3) = obj(:,3)/80 ;

% zhao %maximum 120 wrp
% obj(:,1) = (obj(:,1)/450);
% obj(:,2) = (obj(:,2)/450);
% obj(:,3) = (obj(:,3)/20);
%
%% medieval
% medieval
% obj(:,1) = obj(:,1)*3;
% obj(:,2) = obj(:,2)*3;
% obj(:,3) = obj(:,3)/20;

%%
% three_object
% obj(:,1) = obj(:,1)*3.5;
% obj(:,2) = obj(:,2)*3.5;
% obj(:,3) = obj(:,3)*50;
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

%% OBJECT = BEAR

%{
load bear_67w.txt;
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
obj(:,3) = obj(:,3)/1000;
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
%
% Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
% obj_X = obj(:,1);
% obj_y = obj(:,2);
% obj_z = obj(:,3);
% obj2(:,1) = obj_X(Visible_point_indexes);
% obj2(:,2) = obj_y(Visible_point_indexes);
% obj2(:,3) = obj_z(Visible_point_indexes);
%
% ptCloud = pointCloud(obj2);
% figure,pcshow(ptCloud);
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



%}
%% bear with lesser points
%{
ptCloud = pcread ('../data/bear3.ply');
obj = ptCloud.Location;
obj(:,1)=obj(:,1)+abs(max(obj(:,1))+min(obj(:,1)))/2-0.10;%+abs(min(obj1(:,1)))+1;
obj(:,2)=obj(:,2)+abs(max(obj(:,2))+min(obj(:,2)))/2-0.10;%+abs(min(obj1(:,2)))+1;
% obj(:,3)=obj(:,3)+abs(min(obj(:,3)))/2+1;
obj = obj/300;
ptCloud2 = pointCloud(obj);
ptCloud2.Color=uint8(ptCloud.Color);
figure,pcshow(ptCloud2);
xlabel('X');
ylabel('Y');
zlabel('Z');
obj_c = double(ptCloud2.Color)/256 ;

% obj_c = double(obj_c/256);
%}
%%
max(obj)
min(obj)
% max(obj) - min(obj)
obj_depth = max(obj(:,3)) - min(obj(:,3))
d = 0.7;

file_type = '.bmp';
project_name = mfilename;



Hologram_resolution_x = 1980;
Hologram_resolution_y = 1024;  % Hologram resolution
Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;                     % Hologram resolution
% Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
Hologram_sampling_interval = 7.4e-6; %8e-6;%            % Hologram sampling interval

Nxx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);
Nyy = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );
lambda = [632.8e-9 532e-9 473e-9]; %RGB;


obj_z = obj(:,3);
% obj_z_temp = obj(:,3);
% obj_z= round(obj_z_temp,4);



Cut = sort(unique(obj_z));
% [a,b]=hist(obj_z,unique(obj_z));


ROWS= Hologram_resolution_x;
COLS= Hologram_resolution_y;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));

% total_wrp_list = (105:125);
total_wrp_list = (120);

sub_dir = strcat(project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)));
% mkdir(sub_dir);
% focus_location ='hand_';
% file_name = strcat(focus_location,'time_10_to_',sprintf('%d',Hologram_resolution),'.txt');
% full_file_name = fullfile(sub_dir, file_name);
% fileID = fopen(full_file_name,'w');
prev_wrp_no= -1;


for total_wrp_index = 1 :length(total_wrp_list)
    
    original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    total_wrp = total_wrp_list(total_wrp_index);
    
    size_of_one_depth_range = round(length(Cut) / total_wrp);
    n=numel(Cut);
    m=fix(n/size_of_one_depth_range);
    p=mod(n,size_of_one_depth_range);
    depth_ranges =[mat2cell(Cut(1:m*size_of_one_depth_range),size_of_one_depth_range*ones(m,1),1);{Cut(size_of_one_depth_range*m+1:size_of_one_depth_range*m+p)}];
    
    if(cellfun(@isempty, depth_ranges(end)))
        depth_ranges(end)=[];
    end
    
    if(length(depth_ranges)==prev_wrp_no)
        continue;
    else
        prev_wrp_no = length(depth_ranges);
    end
    
    tic
    for color_index = 1:length(lambda)
        k = 2*pi/lambda(color_index);
        WRPHologram = zeros(Hologram_resolution_x,Hologram_resolution_y);
        counter =0;
        for nwrp = 1:length(depth_ranges)
            WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
            depth_range = cell2mat( depth_ranges(nwrp)) ;
            z_wrp = median(depth_range);  %  WRP location ;
            
            for layer = 1: length(depth_range)
                indexes = find(obj_z==depth_range(layer));
                z =   z_wrp - depth_range(layer);
                N =  (round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 );        %sampling size of N
                Nx = Nxx(indexes);
                Ny = Nyy(indexes);
                color_current_depth_range = obj_c (indexes ,:);
                [y_run, x_run]= meshgrid((-(N-1)/2:(N-1)/2)*Hologram_sampling_interval,(-(N-1)/2:(N-1)/2)*Hologram_sampling_interval);
                r =  ( sign(z)*sqrt(z^2 + y_run.^2 + x_run.^2));
                for n=1:length(indexes)
                    counter = counter +1;
                    fprintf('wrp %d counter: %d color: %d\n',length(depth_ranges), counter ,color_index);
                    %
                    
                    %                       pp = color_current_depth_range(n,color_index)*exp(1i*k*z).*exp(-1i*pi*k*z*r)./r;   %Fourier transform of h
                    %
                    pp = color_current_depth_range(n,color_index)*exp(1j*rand*2*pi)*exp(1j*k*r)./r;
                    
                    WRP(Nx(n):Nx(n)+N-1,Ny(n):Ny(n)+N-1) =  WRP(Nx(n):Nx(n)+N-1,Ny(n):Ny(n)+N-1)+ pp ;
                    
                end
                
            end
            d_wrp_to_hologram = abs(max(d -  z_wrp));
            WRPHologram  = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
            
        end
        
        close all
        for p = 31 %reconstructed
            d2 = d+0.002 - p*0.0010;
            original_red =  FresnelPropogation(k,v, h,-d2,gather(WRPHologram));
            original_abs_red= abs((original_red));
            original(:,:,color_index)=255.*(original_abs_red./max(max(original_abs_red)));
            figure; imshow((abs(original(:,:,color_index)))*10,[]);
            %               toc
            %         file_name = strcat('',project_name,'_color_',num2str(color_index),'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
            %         fullFileName = fullfile(sub_dir,file_name);
            %         save_hologram(gather(WRPHologram),fullFileName);
            
            
        end
        
        %         title(['image o = ',num2str(color_index)])
    end
    
    end_time = toc
    % imshow(rot90(uint8((original))),[]);
    %  uint8(round(RGB64*255));
    % figure;imshow((uint8((original))),[]);
    file_name = strcat(focus_location,project_name,'_d_',num2str(d),'WRP_',num2str(length(depth_ranges)),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    imwrite((uint8((original))),fullFileName);
    fprintf(fileID,' %d %0.3f \n ',length(depth_ranges), end_time);
    
end


%% reconstruction for red
fclose(fileID);