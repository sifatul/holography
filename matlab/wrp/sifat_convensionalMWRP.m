%Main program of hologram generation by wavefront recording plane method
%by Openholo library project
%2017-10-30 update
%

cd D:\Mathlab\wrp
addpath('D:\Mathlab\wrp\data');
clc;clear;close all;
 

load three_obj.txt
obj_name = "three_obj";

obj = zeros(size(three_obj,1),3);
obj(:,1) = three_obj(:,1);
obj(:,2) = three_obj(:,2);
obj(:,3) = three_obj(:,3);

obj(:,1) = (obj(:,1)/max(obj(:,1)))+0.25;
obj(:,2) = obj(:,2)/max(obj(:,2));
obj(:,3) = obj(:,3)/max(obj(:,3));

% obj(:,1) = (obj(:,1)/1000)-min(obj(:,1)/1000);
% obj(:,2) = obj(:,2)/1000-min(obj(:,2)/1000);
% obj(:,3) = obj(:,3)/1000;

obj(:,1) = (obj(:,1)/1000)*2.5; 
obj(:,2) = (obj(:,2)/1000)*2.5;  
obj(:,3) = (obj(:,3)/1000)*20;  
%figure, plot3(obj(:,1),obj(:,2),obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

max(obj) -min(obj)

 obj_c=zeros(length(obj),3);

obj_c(:,1) = three_obj(:,4).*256;
obj_c(:,2) = three_obj(:,5).*256;
obj_c(:,3) = three_obj(:,6).*256; 
ptCloud = pointCloud(obj);
ptCloud.Color=uint8(obj_c);
figure,pcshow(ptCloud); 

%}
%% 
%{

obj_name = 'papilon';

pc = pcread (strcat(obj_name,'.ply'));
pcshow(pc); xlabel('X'); ylabel('Y');zlabel('Z');
obj = im2double(pc.Location);
obj_c = im2double(pc.Color);

 

% temp_idx = find(obj(:,3)<950);
% obj = obj(temp_idx,:);
% ptCloud2 = pointCloud(obj);
% ptCloud2.Color=(pc.Color(temp_idx,:));
% pcshow(ptCloud2);
% % % pcwrite(ptCloud2,'half_medieval');
% 
% obj_2 = im2double(ptCloud2.Location);
% temp_idx2 = find(obj(:,1)<400 );
% obj3 = obj_2(temp_idx2,:);
% ptCloud3 = pointCloud(obj3);
% ptCloud3.Color=(ptCloud2.Color(temp_idx2,:));
% pcshow(ptCloud3);
% obj_c = im2double(ptCloud2.Color);
% figure,pcshow(pc);
% Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
%  
% ptCloud2 = pointCloud(obj(Visible_point_indexes,:));
% ptCloud2.Color=(pc.Color(Visible_point_indexes,:));


% ptCloud2 = pointCloud(obj);
% ptCloud2.Color=(pc.Color);
% pcshow(ptCloud2);xlabel('X');ylabel('Y');zlabel('Z');


obj(:,1) = obj(:,1)/max(obj(:,1)) -0.5;
obj(:,2) = obj(:,2)/max(obj(:,2)) -0.5 ;
obj(:,3) = obj(:,3)/max(obj(:,3));

obj(:,1) = obj(:,1)/1.5;
obj(:,2) = obj(:,2)/1.5;

% obj(:,1) = obj(:,1)/3;
% obj(:,2) = obj(:,2)/3;

obj(:,1) = obj(:,1) /100;
obj(:,2) = obj(:,2) /100;
obj(:,3) = obj(:,3)/80 ;



%}

%% bear with lesser points
 %{
ptCloud = pcread ('bear3.ply');
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
% max(obj)
% min(obj) 
obj_depth = max(obj(:,3)) -min(obj(:,3)) 
d = 0.13;

file_type = '.bmp';
project_name='sifat_conventional_WRP_';


Hologram_resolution = 1024;                       % Hologram resolution    
% Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
Hologram_sampling_interval = 7.4e-6; %8e-6;%            % Hologram sampling interval

Nxx = gpuArray( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution)/2);  
Nyy = gpuArray(round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution)/2 );
lambda = [632.8e-9 532e-9 473e-9]; %RGB; 
 

obj_z = obj(:,3); 
% obj_z_temp = obj(:,3);
% obj_z= round(obj_z_temp,4);


 
Cut = sort(unique(obj_z));
[a,b]=hist(obj_z,unique(obj_z));

ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));

total_wrp_list = [64];
% total_wrp_list = [64];
  
sub_dir = strcat(project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)));
% mkdir(sub_dir);
% file_name = strcat('time_',sprintf('%d',Hologram_resolution),'.txt');
% full_file_name = fullfile(sub_dir, file_name);
% fileID = fopen(full_file_name,'w'); 
 
 
for total_wrp_index = 1 :length(total_wrp_list) 
    
    original=zeros(Hologram_resolution,Hologram_resolution,3); 
    total_wrp = total_wrp_list(total_wrp_index);
    total_wrp
    
    size_of_one_depth_range = round(length(Cut) / total_wrp);
    n=numel(Cut);
    m=fix(n/size_of_one_depth_range);
    p=mod(n,size_of_one_depth_range);
    depth_ranges =[mat2cell(Cut(1:m*size_of_one_depth_range),size_of_one_depth_range*ones(m,1),1);{Cut(size_of_one_depth_range*m+1:size_of_one_depth_range*m+p)}];
    
    if(cellfun(@isempty, depth_ranges(end)))
        depth_ranges(end)=[];
    end
    
     tic
    for color_index = 1:length(lambda)
        k = 2*pi/lambda(color_index); 
        WRPHologram = gpuArray( zeros(Hologram_resolution)); 
        counter =0;
        total = 0; 

        for nwrp = 1:length(depth_ranges) 
            Hologram_wrp = gpuArray( zeros(Hologram_resolution));
            depth_range = cell2mat( depth_ranges(nwrp)) ;
            
            z_wrp = depth_range( round(length(depth_range)/2))+0.00001;  %  WRP location ; 
          

            for layer = 1: length(depth_range)
                indexes = find(obj_z==depth_range(layer)); 
                total = total + 1;               
                z = gpuArray ( abs(z_wrp - depth_range(layer)));                              
                Pobj = obj(indexes,:); 
                N = gpuArray (round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 );        %sampling size of N 
                Nx = Nxx(indexes); 
                Ny = Nyy(indexes);
                color_current_depth_range = obj_c (indexes ,:);
                [y_run, x_run]= meshgrid((-(N-1)/2:(N-1)/2)*Hologram_sampling_interval,(-(N-1)/2:(N-1)/2)*Hologram_sampling_interval);
                r = gpuArray ( sign(z)*sqrt(z^2 + y_run.^2 + x_run.^2));               
               
                for n=1:length(indexes)      
                    counter = counter +1;
                    fprintf('wrp %d counter: %d color: %d\n',total_wrp, counter ,color_index);
                    pp = color_current_depth_range(n,color_index)*exp(1j*2*pi)*exp(1j*k*r)./r;
                    Hologram_wrp(Nx(n):Nx(n)+N-1,Ny(n):Ny(n)+N-1) =  Hologram_wrp(Nx(n):Nx(n)+N-1,Ny(n):Ny(n)+N-1)+ pp ;
                    gather(r);
                end 

            end
            d_wrp_to_hologram = abs(max(d -  z_wrp));            
            WRPHologram  = WRPHologram + FresnelPropagation(gather(Hologram_wrp), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
            
        end

               close all    
        for p = 1:15 %reconstructed   
            d2 = d+0.002 - p*0.0010; 
            original_red =  FresnelPropogation(k,v, h,-d2,WRPHologram);
            original_abs_red= abs(gather(original_red)); 
            original(:,:,color_index)=255.*(original_abs_red./max(max(original_abs_red)));
            figure; imshow((abs(original(:,:,1))),[]);         
        end
        
        title(['image o = ',num2str(color_index)])
    end

end_time = toc
% imshow(rot90(uint8(gather(original))),[]); 
%  uint8(round(RGB64*255));
figure;imshow((uint8(gather(original*2))),[]); 
file_name = strcat('another_focus_full',project_name,'_d_',num2str(d),'threshold_',num2str(total_wrp),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
fullFileName = fullfile(sub_dir,file_name); 
imwrite((uint8(gather(original*2))),fullFileName);
% 
fprintf(fileID,' %d %0.1f \n ',total_wrp, end_time);

end

%% reconstruction for red
  fclose(fileID);


toc
% file_prefix = 'hologram_red';
% file_name = strcat(obj_name,project_name,sprintf('%0.1f',Hologram_resolution),sprintf('%0.9f',Hologram_sampling_interval),file_prefix,file_type);
% save_hologram(WRPHologram_red,file_name);
imshow(rot90(uint8(gather(original))),[]);