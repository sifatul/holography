addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
clc;clear;close all;
papilon
obj(:,1) = (obj(:,1)-0.0005)/2.5;
obj(:,2) = (obj(:,2)-0.0006)/2.5;
obj(:,3) = obj(:,3)/18 ;

obj_depth = max(obj) - min(obj)


file_type = '.bmp';
project_name = mfilename;


Hologram_resolution_x = 1280;
Hologram_resolution_y = 720;  % Hologram resolution
Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;
Hologram_sampling_interval = 8e-6; %8.4e-6;%7.4e-6;8.4e-6; %;            % Hologram sampling interval



% Hologram_wrp  =( zeros(Hologram_resolution, Hologram_resolution, 3));  %( zeros(Hologram_resolution));

lambda = [632.8e-9 532e-9 473e-9]; %RGB;

k = 2*pi./lambda;

obj_z = obj(:,3);

Cut = sort(unique(obj_z));

ROWS= Hologram_resolution_x;
COLS= Hologram_resolution_y;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));


obj_depth= max(obj(:,3)) -min(obj(:,3));
d = 0.25;

Nx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);
Ny = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );

max_counter= 0;


sub_dir = strcat('D:\Mathlab\wrp\final',obj_name,'',sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)),'\',project_name);
mkdir(sub_dir );
%
file_name = strcat('new_time_temp_',sprintf('%d',Hologram_resolution),'.txt');
full_file_name = fullfile(sub_dir, file_name);
fileID = fopen(full_file_name,'w');

for ii=1:1
    
    
    original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    original2=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
          WRPHologram =   zeros(Hologram_resolution_x,Hologram_resolution_y,3);
  
    tic
    for color_index = 1:length(lambda)
        %     color_index
        counter = 0;
        WRP =  (zeros(Hologram_resolution_x,Hologram_resolution_y));
  
        counter2 =0;
        
        for o = 1: length(Cut)
            indexes = find(obj_z==Cut(o));
            
            z_wrp = obj_z(indexes)+0.005;  %  WRP location
            WRP =  zeros(Hologram_resolution_x,Hologram_resolution_y);
           
            z = (z_wrp -  Cut(o));
            color_for_this_depth = obj_c(indexes,color_index);
            
            %         active_area =(round(abs(lambda(color_index).*z(1)./(Hologram_sampling_interval^2)/2)+0.5).*2-1 );        %sampling size of N
            active_area = round(2*abs(z).*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
            if (min(active_area) <= 0)
                return;
            end
            Nx_for_this_depth = Nx(indexes);
            Ny_for_this_depth = Ny(indexes);
            
            xn= active_area;yn = active_area;
            nxn=Nx_for_this_depth+active_area; nyn=Ny_for_this_depth+active_area;
            if (nxn>Hologram_resolution_x)
                xn=Hologram_resolution_x-Nx_for_this_depth;
                nxn=Hologram_resolution_x;
            end
            if(nyn>Hologram_resolution_y)
                yn=Hologram_resolution_y-Ny_for_this_depth;
                nyn=Hologram_resolution_y;
            end
            
            
            [y_run, x_run]= meshgrid((-(active_area-1)/2:(active_area-1)/2)*Hologram_sampling_interval,(-(active_area-1)/2:(active_area-1)/2)*Hologram_sampling_interval);
            r =  (sign(z(1))*sqrt(z(1)^2 + y_run.^2 + x_run.^2));
            for Ni = 1 : length( indexes)
                nx=Nx_for_this_depth(Ni); ny=Ny_for_this_depth(Ni);
                counter = counter +1;
                fprintf('counter: %d color: %d\n',counter ,color_index);
                Sub_hologram =  color_for_this_depth(Ni)*exp(-1j*2*rand*pi)*exp(1j*k(color_index)*r)./r ;
                WRP(Nx_for_this_depth:nxn-1,Ny_for_this_depth:nyn-1) =  WRP(Nx_for_this_depth:nxn-1,Ny_for_this_depth:nyn-1)+ Sub_hologram;
                
            end
            d_wrp_to_hologram = max(abs(d-z_wrp));
            WRPHologram(:,:,color_index) = WRPHologram(:,:,color_index) + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
            
            
        end %end of all layers
        
        
        
        
        %   disp( color_index);
        
        %     title(num2str(color_index));
    end
    
    close all
    for p = 17%reconstructed
        for c=1:3
            d2 = d+0.002 - p*0.001
            original_red =  FresnelPropogation(k(c),v, h,-d2,WRPHologram(:,:,c));
            original_abs_red= abs((original_red));
            original(:,:,c)=original_abs_red;
            original2(:,:,c)= 255.*(original_abs_red./max(max(original_abs_red))) ;
            figure; imshow(rot90(original_abs_red),[]);
            
            %
            %      file_name = strcat('table_box_hologram',project_name,'_color_',num2str(color_index),'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
            %     fullFileName = fullfile(sub_dir,file_name);
            %     save_hologram((WRPHologram),fullFileName);
        end
         figure; imshow(rot90(original),[]);
    end
    
    %% reconstruction for red
    end_time= toc;
    imshow(rot90(uint8((original2*2 ))),[]);
    
    
    % fprintf(fileID,' %d %0.1f \n ',total_wrp, end_time);
    % file_name = strcat('full_hologram',project_name,'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    % fullFileName = fullfile(sub_dir,file_name);
    % save_hologram((WRPHologram),fullFileName);
    % % imshow(rot90(uint8((original))),[]);
    file_name = strcat('book2',project_name,'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    imwrite((uint8((original2*2))),fullFileName);
    
    fprintf(fileID,' %d %0.1f \n ',total_wrp, end_time);
    %  disp(counter2_limit);
end

fclose(fileID);
