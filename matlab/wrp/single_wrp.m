addpath('D:\Mathlab\wrp\data');
addpath('D:\My Research\holography\matlab\wrp\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
clc;clear;close all;

% me_cube2 
% obj(:,1) = ((obj(:,1) )/400);
% obj(:,2) = ((obj(:,2)+0.15 )/350);
% obj(:,3) = (obj(:,3))/40;
% 
% max(obj)
% min(obj)
% max(obj) - min(obj)
% obj_depth = max(obj(:,3)) - min(obj(:,3))
% d = 0.5;


% zhao %maximum 120 wrp
% obj(:,1) = (obj(:,1)/450);
% obj(:,2) = (obj(:,2)/450);
% obj(:,3) = (obj(:,3)/20);

four_items
% obj(:,1) = ((obj(:,1) )/600);
% obj(:,2) = ((obj(:,2)  )/1000);
% obj(:,3) = (obj(:,3))/30;
% max(obj)
% min(obj)
% 
%1980*1024
obj(:,1) = ((obj(:,1) )/300);
obj(:,2) = ((obj(:,2)  )/500);
obj(:,3) = (obj(:,3))/30;
max(obj)
min(obj)



obj_depth = max(obj(:,3)) - min(obj(:,3))
d = 0.25 ;
file_type = '.bmp';
project_name = mfilename;



Hologram_resolution_x = 1980;
Hologram_resolution_y = 1024;  % Hologram resolution
Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;
Hologram_sampling_interval = 7.4e-6; %8e-6;%            % Hologram sampling interval

Nxx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);
Nyy = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );
lambda = [632.8e-9 532e-9 473e-9]; %RGB;


obj_z = obj(:,3);
Cut = sort(unique(obj_z));


ROWS= Hologram_resolution_x;
COLS= Hologram_resolution_y;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));

% total_wrp_list = (105:125);

counter = 0;
original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
for color_index = 1:length(lambda)  %iterate each color
    k = 2*pi/lambda(color_index);
    WRPHologram =  zeros(Hologram_resolution_x,Hologram_resolution_y);
    
    for point = 1: length(obj)
        
        WRPHologram  = WRPHologram + FresnelPropagation(obj(point,:), Hologram_sampling_interval, Hologram_sampling_interval, d, lambda(color_index));
        
         
      
    end
    
    
    
    
    close all
    for p = 22   %reconstructed
        d2 = d+0.002 - p*0.0010
        original_red =  FresnelPropogation(k,v, h,-d2,gather(WRPHologram));
        original_abs_red= abs((original_red));
        original(:,:,color_index)=255.*(original_abs_red./max(max(original_abs_red)));
                     figure; imshow(rot90(abs(original(:,:,color_index)))*10,[]);title(p);
        %                           toc
        %             file_name = strcat(project_name,'_color_',num2str(color_index),'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
        %             fullFileName = fullfile(sub_dir,file_name);
        %                     save_hologram(gather(WRPHologram),fullFileName);
        %              imwrite((uint8((original))),fullFileName);
        %           imwrite(uint8(rot90((original(:,:,color_index)))*10),fullFileName);
        
    end
end

 imshow(rot90(uint8((original*2))),[]);
sub_dir = strcat(project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)));
mkdir(sub_dir);

    file_name = strcat('boook',project_name,'_d_',num2str(d), '_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    imwrite((uint8((original*2))),fullFileName);
