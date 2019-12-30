addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
cd 'D:\My Research\holography\matlab\wrp'
clc;clear;close all;clear all;

four_items
% four_items_new
obj(:,1) = ((obj(:,1) )/800); %%800
obj(:,2) = ((obj(:,2)  )/800);
obj(:,3) = (obj(:,3))/10;




%%
max(obj)
min(obj)
% max(obj) - min(obj)
obj_depth = max(obj(:,3)) - min(obj(:,3))
% d = 0.18;
d = 0.15;
file_type = '.bmp';

%% Hologram parameter

Hologram_resolution_x = 720;
Hologram_resolution_y = 720;  % Hologram resolution
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
[Cut,~,idx] = unique(obj_z);
n = accumarray(idx(:),1);


his =  histogram(Cut,19)
BinCounts = his.BinCounts;
Values =his.Values;
BinEdges = his.BinEdges;
start = 1;
tail =1;
t = 1;
C = {};
for i = 1:length(BinCounts)
    tail=i;
    if BinCounts(i)== 0
        tail=i; 
        if  start~= tail 
            temp  = BinEdges(start:tail-1);
            %%depth ranges based on bin edges unfinished
            C{t}= temp;            
            start = i; 
        end
        
    end 
end
if start~= tail
    C{t}= temp;
end

% WRP_list= [1,10,20,30,50,90,120,189];
WRP_list= [8];
% WRP_list= [2;4;8;16;32;66];
prev_wrp_no = -1;

for i= 1: length(WRP_list)
    
    original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    total_wrp = WRP_list(i);
    
    % % % %     %     conventional M-WRP Depth range segmentation
    %     project_name = mfilename;
    %     depth_ranges = depth_segmentation_conventional(Cut,total_wrp);
    %
    
    % load DM-WRP segmentation
    core_dmwrp;
    %
    %% data storage
    
    if(length(depth_ranges)==prev_wrp_no)
        continue;
    else
        prev_wrp_no = length(depth_ranges);
    end
    
    if length(cell2mat(depth_ranges(1))) ~= length(cell2mat(depth_ranges(end)))
        disp("#Unequal depth ranges#");
    end
    
    
    
    %
    %     sub_dir = strcat('D:\Mathlab\wrp\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)),'\',project_name);
    %     focus_location ='cactus';
    %     mkdir(sub_dir);
    %     file_name = strcat('time','.txt');
    %     full_file_name = fullfile(sub_dir, file_name);
    %     fileID = fopen(full_file_name,'a');
    
    
    
    %
    
    tic
    for color_index = 1:length(lambda)
        layer = 1;
        k = 2*pi/lambda(color_index);
        WRPHologram = zeros(Hologram_resolution_x,Hologram_resolution_y);
        counter =1;
        for nwrp = 1:length(depth_ranges)
            WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
            depth_range = cell2mat( depth_ranges(nwrp)) ;
            %             if nwrp == length(depth_ranges)
            %                 depth_ranges(end) =[];
            %             end
            min_z = min(depth_range);
            max_z = max(depth_range);
            
            if mod(length(depth_range),2) == 1
                depth_range(end+1)=depth_range(end)+0.005;
                z_wrp = median(depth_range);
                depth_range(end) = [];
            else
                z_wrp = median(depth_range);  %  WRP location ;
            end
            %             if length(depth_range)== 1
            %                 z_wrp = z_wrp + 0.000005;
            %             elseif find(Cut == z_wrp) > 0
            %                 z_wrp = z_wrp + min(unique( diff(depth_range)))/2;
            %             end
            %
            indexes = find(obj_z>= min_z & obj_z <= max_z);
            
            z = abs(z_wrp-obj_z(indexes));
            %             N = round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 ;        %sampling size of N
            %
            N = round(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
            %                         for p = 1 :length(depth_range)
            %                             zp = abs(z_wrp-depth_range(p));
            %                             active(layer) = round(zp.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval) ;
            %                             layer = layer+1;
            %                         end
            % % %
            %             if min(N) ==0
            %                 return
            %             end
            Nx = Nxx(indexes);
            Ny = Nyy(indexes);
            color_current_depth_range = obj_c(indexes ,:);
            
            
            for index = 1: length(indexes)
                
                [y_run, x_run]= meshgrid((-(N(index)-1)/2:(N(index)-1)/2)*Hologram_sampling_interval,(-(N(index)-1)/2:(N(index)-1)/2)*Hologram_sampling_interval);
                r = sign(z(index))*sqrt(z(index)^2 + y_run.^2 + x_run.^2);
                pp = color_current_depth_range(index,color_index)*exp(-1j*2*rand*pi)*exp(1j*k*r)./r ;
                
                WRP(Nx(index):Nx(index)+N(index)-1,Ny(index):Ny(index)+N(index)-1) =  WRP(Nx(index):Nx(index)+N(index)-1,Ny(index):Ny(index)+N(index)-1)+ pp;
                counter = counter + 1 ;
                fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, length(depth_ranges)); %length(depth_ranges)
                
            end
            
            d_wrp_to_hologram = abs(d -  z_wrp);
            WRPHologram  = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
            
        end
        
        
        close all
        for p =  46 %reconstructed
            d2 = d+0.002 - p*0.0010
            original_red =  FresnelPropogation(k,v, h,-d2,gather(WRPHologram));
            original_abs_red= abs((original_red));
            original(:,:,color_index)=255.*(original_abs_red./max(max(original_abs_red)));
            %                         figure; imshow(rot90(abs(original(:,:,color_index)))*10,[]);title(p);
            %             toc
            %             file_name = strcat(project_name,'_color_',num2str(color_index),'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
            %             fullFileName = fullfile(sub_dir,file_name);
            %                     save_hologram(gather(WRPHologram),fullFileName);
            %              imwrite((uint8((original))),fullFileName);
            %           imwrite(uint8(rot90((original(:,:,color_index)))*10),fullFileName);
            
        end
        
        %
    end
    
    end_time = toc;
    fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, length(depth_ranges)); %length(depth_ranges)
    % imshow(rot90(uint8((original*2))),[]);
    %  uint8(round(RGB64*255));
    % figure;imshow((uint8((original))),[]);
    file_name = strcat('cup2','_d_',num2str(d),'WRP_',num2str(length(depth_ranges)),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    imwrite((uint8((original*2))),fullFileName);
    fprintf(fileID,' %d %0.3f \n ',length(depth_ranges), end_time);
    
end


%% reconstruction for red
fclose(fileID);
x = time(:,1);
y = time(:,2);
plot( x, y);
idx = islocalmin(y);
figure(1)
hold on
plot(x,y)
plot(x(idx),y(idx),'*r')
legend('Curve','Local Min')
hold off
legend({'cos(x)','cos(2x)','cos(3x)','cos(4x)'},'Location','northwest','NumColumns',2)
legend('Conventional M-WRP','Proposed DM-WRP','Location','northeast')
% title('Segmentation in Conventional MWRP');
ylabel('distance in Z-axis(m)')
xlabel('No. of depth layers')
hold on;
