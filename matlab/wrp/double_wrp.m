addpath('D:\Mathlab\wrp\data');
addpath('D:\My Research\holography\matlab\wrp\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
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
three_object  
obj(:,1) = obj(:,1)*3.5; 
obj(:,2) = obj(:,2)*3.5;  
obj(:,3) = obj(:,3)*50;   

max(obj)
min(obj)
% max(obj) - min(obj)
obj_depth = max(obj(:,3)) - min(obj(:,3))
d = 0.15;

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
total_wrp_list = (2:105);
% saving_dir = 'D:\Mathlab\wrp\';
% sub_dir = strcat(project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)));
% mkdir(sub_dir);
% focus_location ='hand_';
% file_name = strcat(focus_location,'time_10_to_',sprintf('%d',Hologram_resolution),'.txt');
% full_file_name = fullfile(sub_dir, file_name);
% fileID = fopen(full_file_name,'w'); 
prev_wrp_no= -1;
   original=zeros(Hologram_resolution_x,Hologram_resolution_y,3); 
 for color_index = 1:length(lambda)  %iterate each color
    k = 2*pi/lambda(color_index);   
    WRPHologram =  zeros(Hologram_resolution_x,Hologram_resolution_y);
    WRPHologram2 =  zeros(Hologram_resolution_x,Hologram_resolution_y);
    
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
    
 end


