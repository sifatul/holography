addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
cd 'D:\My Research\holography\matlab\wrp'
clc;clear;close all;clear all;
% papilon
% obj(:,1) = (obj(:,1)-0.0005)/2.5;
% obj(:,2) = (obj(:,2)-0.0006)/2.5;
% obj(:,3) = obj(:,3)/10 ;

% four_items_new
% obj(:,1) = ((obj(:,1)-0.1 )/150);
% obj(:,2) = ((obj(:,2)+1.65  )/1800);
% obj(:,3) = (obj(:,3))/6.5;
% three_object
% %
% obj(:,1) = (obj(:,1)+0.0002).*4  ;
% obj(:,2) = (obj(:,2)-0.0002)*4  ;
% obj(:,3) = obj(:,3)*150;


% test_three_item
% obj(:,1) = ((obj(:,1)+0.04 )/200);
% obj(:,2) = ((obj(:,2)+1.8  )/1500);
% obj(:,3) = (obj(:,3))/6.5;

%
%%
three_object

  

for it=1:3
    
    
    if (it==1)
        component_idx = find(  obj_temp(:,3)<0.0007);
    elseif (it==2)
        component_idx = find( obj_temp(:,3)>=0.0007 & obj_temp(:,3)<=0.0009);
    else
        component_idx = find(  obj_temp(:,3)>=0.0009);
    end
    
    obj = obj_temp(component_idx,:);
    obj_c = obj_temp_c(component_idx,:);
    obj(:,1) = (obj(:,1)+0.0003)*2.5 ;
    obj(:,2) = obj(:,2) *2.5;
    obj(:,3) = obj(:,3)*110;
    
    
    
    
    
    max(obj)
    min(obj)
    %
    % % max(obj) - min(obj)
    obj_depth = max(obj_temp(:,3)) - min(obj_temp(:,3))
    % d = 0.35;
    d = 0.35;
    
    Hologram_resolution_x = 1024;
    Hologram_resolution_y = 1024  % Hologram resolution
    
    Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;                     % Hologram resolution
    % Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
    Hologram_sampling_interval = 8e-6; %8e-6;%            % Hologram sampling interval
    
    Nxx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);
    Nyy = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );
    
    
    
    lambda = [632.8e-9 532e-9 473e-9]; %RGB;
    k = 2*pi./lambda;
    
    %% Fresnel propagation field
    
    ROWS= Hologram_resolution_x;
    COLS= Hologram_resolution_y;
    v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
    h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
    
    
    %%
    obj_z = (obj(:,3));
    Cut = sort(unique(obj_z));
    layer_diff = abs(Cut(1)-Cut(2));
    [~,~,idx] = unique(obj_z);
    obj_point = accumarray(idx(:),1);
    
    figure;scatter([1:length(Cut)],Cut);
    % figure;scatter([1:length(range)],range);
 
    
    project_name = mfilename;
    
    sub_dir = strcat('D:\sifatul\',obj_name,'',sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(unique(obj_temp(:,3)))),'\',project_name);
    
    focus_location ='fly';
    mkdir(sub_dir);
    file_name = strcat('time_cube',num2str(it),'_', datestr(now, 'dd-HH-MM'),'.txt');
    time_file = fullfile(sub_dir, file_name);
    fileID = fopen(time_file,'a');
    prev_wrp_no= -1;
    
    WRP_list = [1:20];
    WRPHologram = zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    for i= 1: length(WRP_list)
        
        
        original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
        
        active_area_limit =  WRP_list(i);
        
        for color_index = 1:length(lambda)
            depth_ranges = [];
            %         active_area_limit =  WRP_list(color_index);
            ztan = (active_area_limit.*Hologram_sampling_interval)./(2*tan(lambda(color_index) ./(2.*Hologram_sampling_interval)));
            
            z_wrp = (ztan+Cut(1));
            wrp_no = 1;
            temp_range = [];
            C1 ={};
            t = 0;
            for layer = 1: length(Cut)
                z = abs(z_wrp-Cut(layer));
                N = round(2*abs(z).*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
                
                if  N <= active_area_limit
                    t = t+1;
                    temp_range(t) = Cut(layer);
                elseif N > active_area_limit
                    C1{wrp_no} = temp_range;
                    temp_range = [];
                    t = 1;
                    wrp_no = wrp_no + 1 ;
                    z_wrp = abs(ztan+Cut(layer));
                    temp_range(1) = Cut(layer);
                    %             z = abs(z_wrp-Cut(layer));
                    %             N = round(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
                end
            end
            C1{wrp_no} = temp_range;
            
            depth_ranges = C1;
            
            tic
            
            
            counter =1;
            for nwrp = 1:length(depth_ranges)
                WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
                depth_range = cell2mat( depth_ranges(nwrp)) ;
                %             min_z = min(depth_range);
                %             max_z = max(depth_range);
                
                
                if mod(length(depth_range),2) == 1
                    depth_range(end+1)=depth_range(end)+layer_diff;
                    z_wrp = median(depth_range);
                    depth_range(end) = [];
                else
                    z_wrp = median(depth_range);  %  WRP location ;
                end
                
                for layer = 1: length(depth_range)
                    indexes = find(obj_z== depth_range(layer));
                    
                    z = max((z_wrp-obj_z(indexes)));
                    %             N = round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 ;        %sampling size of N
                    %
                    N = round(2*abs(z).*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
                    
                    %                                 if N==0
                    %                                     return;
                    %                                 end
                    Nx = Nxx(indexes);
                    Ny = Nyy(indexes);
                    
                    
                    color_current_depth_range = obj_c(indexes ,color_index);
                    [y_run, x_run]= meshgrid((-(N-1)/2:(N-1)/2)*Hologram_sampling_interval,(-(N-1)/2:(N-1)/2)*Hologram_sampling_interval);
                    r = sign(z)*sqrt(z^2 + y_run.^2 + x_run.^2);
                    
                    
                    for index = 1: length(indexes)
                        xn= N;yn = N;
                        nxn=Nx(index)+N; nyn=Ny(index)+N;
                        if (nxn>Hologram_resolution_x)
                            xn=Hologram_resolution_x-Nx(index);
                            nxn=Hologram_resolution_x;
                        end
                        if(nyn>Hologram_resolution_y)
                            yn=Hologram_resolution_y-Ny(index);
                            nyn=Hologram_resolution_y;
                        end
                        
                        Sub_hologram = color_current_depth_range(index)*exp(1j*rand*2*pi)*exp(1j*  k(color_index)*r(1:xn,1:yn))./r(1:xn,1:yn) ;
                        WRP(Nx(index):nxn-1,Ny(index):nyn-1) =  WRP(Nx(index):nxn-1,Ny(index):nyn-1)+ Sub_hologram;
                        
                        %                     Sub_hologram = color_current_depth_range(index,color_index)*exp(-1j*2*rand*pi)*dir ;
                        
                        
                        counter = counter + 1 ;
                        fprintf('counter: %d color: %d active_area_limit: %d \n',counter ,color_index, active_area_limit); %length(depth_ranges)
                        
                    end
                    
                    
                end
                d_wrp_to_hologram = abs(d -  z_wrp);
                WRPHologram(:,:,color_index)  = WRPHologram(:,:,color_index) + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
                
            end
            fprintf(fileID,' %d %d %d %0.3f \n ',length(depth_ranges),active_area_limit, color_index, toc);
            %         end_time(color_index) = toc;
            
        end
        
        
        
    end
    
    
    
    fclose(fileID);
    
end

sub_dir

sub_dir





% c = 1;
% for i= 1:length(R)
%     if mod( R(i,2),2)==0
%         time2(c,:)= R(i,:);
%         c = c+1;
%     end
% end
imid_data = load(time_file);
for t=1:3
    R = imid_data(find(imid_data(:,3)== t ),:);
    % time = imid(find( mod(R(:,2),2)==0),:)
    [val,idx] = min(R(:,4));
    
    output(t,:) =  R(idx,:);
    %     R(idx,:)
    time(t) = R(idx,end);
end
output
time'
sum(time)