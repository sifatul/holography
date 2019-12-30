addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
cd 'D:\My Research\holography\matlab\wrp'
clc;clear;close all;clear all;

four_items_new
obj(:,1) = ((obj(:,1) )/800);
obj(:,2) = ((obj(:,2)  )/800);
obj(:,3) = (obj(:,3))/10;


%%
max(obj)
min(obj)

% max(obj) - min(obj)
obj_depth = max(obj(:,3)) - min(obj(:,3))
d = 0.18;
file_type = '.bmp';

%% Hologram parameter

Hologram_resolution_x = 720;
Hologram_resolution_y = 720;
Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;                     % Hologram resolution
% Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
Hologram_sampling_interval = 7.4e-6; %8e-6;%            % Hologram sampling interval

Nxx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);
Nyy = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );



lambda = [632.8e-9 532e-9 473e-9]; %RGB;

%% Fresnel propagation field

ROWS= Hologram_resolution_x;
COLS= Hologram_resolution_y;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));


%%


obj_z = (obj(:,3));
Cut = sort(unique(obj_z));

active_area_limit =  5;
ztan = (active_area_limit.*Hologram_sampling_interval)./tan(lambda./(2.*Hologram_sampling_interval));
for color_index = 1:length(lambda)
    z_r = ztan(color_index);
    z_wrp = (z_r+Cut(1));
    wrp_no = 1;
    temp_range = {};
    C1 ={};
    t = 0;
    for layer = 1: length(Cut)
        z = abs(z_wrp-Cut(layer));
        N = round(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
        
        if  N <= active_area_limit
            t = t+1;
            temp_range{t} = Cut(layer);
        elseif N > active_area_limit
            C1{wrp_no} = temp_range;
            temp_range = {};
              t = 1;
              wrp_no = wrp_no + 1 ;
            z_wrp = abs(z_r+Cut(layer));
             temp_range{1} = Cut(layer);
%             z = abs(z_wrp-Cut(layer));
%             N = round(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
        end
    end
    
    toc;
    
    
end









% project_name = mfilename;
% sub_dir = strcat('D:\Mathlab\wrp\',project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)));
% focus_location ='cup';
% mkdir(sub_dir);
% file_name = strcat('time_hope','.txt');
% full_file_name = fullfile(sub_dir, file_name);
% fileID = fopen(full_file_name,'a');
% prev_wrp_no= -1;

WRP_list = [1:20];
% depth_ranges = zeros(3);

for i= 1: length(WRP_list)
    
    original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    
    %%%% conventional M-WRP Depth range segmentation
    project_name = mfilename;
    active_area_limit =  WRP_list(i);
    ztan = (active_area_limit.*Hologram_sampling_interval)./tan(lambda./(2.*Hologram_sampling_interval));
    tic
    for color_index = 1:length(lambda)
        k = 2*pi/lambda(color_index);
        WRPHologram = zeros(Hologram_resolution_x,Hologram_resolution_y);
        WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
        counter =0;
        z_r = ztan(color_index);
        z_wrp = (z_r+Cut(1));
%         WRP_counter =1;
        
        for layer = 1: length(Cut)
            indexes = find(obj_z == Cut(layer));
            z = abs(z_wrp-Cut(layer));
            N = round(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
            
            if N > active_area_limit  || N == 0
%                 WRP_counter = WRP_counter+1;
                d_wrp_to_hologram = abs(max(d -  z_wrp));
                WRPHologram  = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
                
                WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
                z_wrp = abs(z_r+Cut(layer));
                z = abs(z_wrp-Cut(layer));
                N = round(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
            end
            
            Nx = Nxx(indexes);
            Ny = Nyy(indexes);
            color_current_depth_range = obj_c(indexes ,:);
            [y_run, x_run]= meshgrid((-(N-1)/2:(N-1)/2)*Hologram_sampling_interval,(-(N-1)/2:(N-1)/2)*Hologram_sampling_interval);
            r = sign(z)*sqrt(z^2 + y_run.^2 + x_run.^2);
            bb = exp(1j*k*r)./r;
            for index = 1: length(indexes)
                
                pp = color_current_depth_range(index,color_index)* exp(-1j*2*rand*pi)*bb;
                WRP(Nx(index):Nx(index)+N-1,Ny(index):Ny(index)+N-1) =  WRP(Nx(index):Nx(index)+N-1,Ny(index):Ny(index)+N-1)+ pp;
                counter = counter + 1 ;
%                                 fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, N); %length(depth_ranges)
                
            end
            
            
        end
        
        d_wrp_to_hologram = abs(max(d -  z_wrp));
        WRPHologram  = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
        
        
        close all
        for p = 94 %reconstructed
            d2 = d+0.002 - p*0.0010
            original_red =  FresnelPropogation(k,v, h,-d2,gather(WRPHologram));
            original_abs_red= abs((original_red));
            original(:,:,color_index)=255.*(original_abs_red./max(max(original_abs_red)));
%                         figure; imshow(rot90(abs(original(:,:,color_index)))*10,[]);title(p);
            %                           toc
            %             file_name = strcat(project_name,'_color_',num2str(color_index),'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
            %             fullFileName = fullfile(sub_dir,file_name);
            %                     save_hologram(gather(WRPHologram),fullFileName);
            %              imwrite((uint8((original))),fullFileName);
            %                       imwrite(uint8(rot90((original(:,:,color_index)))*10),fullFileName);
            
        end
        
        end_time(color_index) = toc;
        % imshow(rot90(uint8((original*2))),[]);
        
    end
    
    
    % imshow(rot90(uint8((original*2))),[]);
    %  uint8(round(RGB64*255));
    % figure;imshow((uint8((original))),[]);
    fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, N); %length(depth_ranges)
    file_name = strcat('cup_hope','_d_',num2str(d),'_active_area_limit_',num2str(N),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    imwrite((uint8((original*2))),fullFileName);
    
    
    fprintf(fileID,' %d %0.3f \n ',N, sum(end_time));
    %     end_time
end


fclose(fileID);