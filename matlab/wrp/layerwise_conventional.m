addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
cd 'D:\My Research\holography\matlab\wrp'
clc;clear;close all;   clear all;

% papilon
% obj(:,1) = obj(:,1)/2.5;
% obj(:,2) = obj(:,2)/2.5;
% obj(:,3) = obj(:,3)/20 ;

% zhao %maximum 120 wrp
% zhao
% obj(:,1) = (obj(:,1)/450);
% obj(:,2) = (obj(:,2)/450);
% obj(:,3) = (obj(:,3)/20);
%% medieval
% medieval
% obj(:,1) = obj(:,1)*3;
% obj(:,2) = obj(:,2)*3;
% obj(:,3) = obj(:,3)/20;

%%
% three_object
% % obj(:,1) = obj(:,1)*3.5;
% % obj(:,2) = obj(:,2)*3.5;
%
% obj(:,1) = obj(:,1);
% obj(:,2) = obj(:,2);
% obj(:,3) = obj(:,3)*90;
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

%  four_items
% % % % 65 segmentation
% obj(:,1) = ((obj(:,1) )/600);
% obj(:,2) = ((obj(:,2)  )/1000);
% obj(:,3) = (obj(:,3))/30;
% max(obj)
% min(obj)
%
%
% % %
% % 1980*1024
% obj(:,1) = ((obj(:,1) )/300);
% obj(:,2) = ((obj(:,2)  )/500);
% obj(:,3) = (obj(:,3))/30;
% max(obj)
% min(obj)



% one_book_cube
% obj(:,1) = ((obj(:,1) )/500);
% obj(:,2) = ((obj(:,2)+0.25 )/600);
% obj(:,3) = (obj(:,3))/40;
%
%
% max(obj)
% min(obj)
% me_cube
% obj(:,1) = ((obj(:,1) )/450);
% obj(:,2) = ((obj(:,2)+0.15 )/600);
% obj(:,3) = (obj(:,3))/20;
%
% %
% book_cube
% obj(:,1) = ((obj(:,1) )/150);
% obj(:,2) = ((obj(:,2)+0.45 )/2000);
% obj(:,3) = (obj(:,3))/20;




me_cube2
obj(:,1) = ((obj(:,1) +0.1)/520);
obj(:,2) = ((obj(:,2)+0.15 )/670);
obj(:,3) = (obj(:,3))/22;
% %

max(obj)
min(obj)


%1980*1024
% me_cube2
% obj(:,1) = ((obj(:,1) )/400);
% obj(:,2) = ((obj(:,2)+0.15 )/350);
% obj(:,3) = (obj(:,3))/40;


% max(obj)
% min(obj)

% cup_cube
% obj(:,1) = ((obj(:,1) )/750);
% obj(:,2) = ((obj(:,2)+0.50)/900);
% obj(:,3) = (obj(:,3))/20;
% max(obj)
% min(obj)

% me_cube
% obj(:,1) = ((obj(:,1) )/450);
% obj(:,2) = ((obj(:,2)+0.15 )/600);
% obj(:,3) = (obj(:,3))/40;
%
%
% max(obj)
% min(obj)

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
 obj(:,1) =  obj(:,1)./10 ;
  obj(:,2) =  obj(:,2)./10 ;

Visible_point_indexes = HPR(obj, max(obj(:,2)),5);
obj_X = obj(:,1);
obj_y = obj(:,2);
obj_z = obj(:,3);
obj2(:,1) = obj_X(Visible_point_indexes);
obj2(:,2) = obj_y(Visible_point_indexes);
obj2(:,3) = obj_z(Visible_point_indexes);

ptCloud = pointCloud(obj2);
ptCloud.Color = obj_c(Visible_point_indexes,:);
figure,pcshow(ptCloud);
 



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
Cut = sort(unique(obj_z));
% obj_z_temp = obj(:,3);
% obj_z= round(obj_z_temp,4);



% Cut = sort(unique(obj_z));

% [a,b]=hist(obj_z,unique(obj_z));

% N = length(Cut);
% x = 1:N;
% WRP_list = x(~(rem(N, x)));
% WRP_list(end) = 248;
% WRP_list(1) = [];
WRP_list= [9];
prev_wrp_no= -1;

for i= 1: length(WRP_list)
    
    original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    
    %% % % conventional M-WRP Depth range segmentation
    project_name = mfilename;
    total_wrp = WRP_list(i);
    depth_ranges = depth_segmentation_conventional(Cut,total_wrp);
    if(length(depth_ranges)==prev_wrp_no)
        continue;
    else
        prev_wrp_no = length(depth_ranges);
    end
    
    %% load DM-WRP segmentation
    %     core_dmwrp;
    
    %% data storage
    % t = 0;
    % for nwrp = 1:length(depth_ranges)
    %         depth_range = cell2mat( depth_ranges(nwrp)) ;
    %         if length( diff (depth_range))>= 1
    %             t = t + sum(diff (depth_range))
    %         end
    %
    % end
    % t/nwrp
    if length(cell2mat(depth_ranges(1))) ~= length(cell2mat(depth_ranges(end)))
        disp("#Unequal depth ranges#");
    end
    
    
    
    %     sub_dir = strcat(project_name,'\',obj_name,sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)),'\final');
    %     focus_location ='person';
    %     mkdir(sub_dir);
    %     file_name = strcat('time_12','.txt');
    %     full_file_name = fullfile(sub_dir, file_name);
    %     fileID = fopen(full_file_name,'a');
    %     prev_wrp_no= -1;
    % %
    
    
    tic
    for color_index = 1:length(lambda)
        k = 2*pi/lambda(color_index);
        WRPHologram = zeros(Hologram_resolution_x,Hologram_resolution_y);
        counter =1;
        t=1;
        for nwrp = 1:length(depth_ranges)
            WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
            depth_range = cell2mat( depth_ranges(nwrp)) ;
            z_wrp = median(depth_range);  %  WRP location ;
            if find(Cut == z_wrp) > 0
                z_wrp = z_wrp + min(unique( diff(depth_range)))/2;
            end
%             total =0;
            for layer = 1: length(depth_range)
                indexes= find(obj_z== depth_range(layer));
%                 total = total + length(indexes);
                z = unique(abs(z_wrp-obj_z(indexes)));
                N = round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 ;        %sampling size of N
                %             N  = ceil(z.*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
                Nx = Nxx(indexes);
                Ny = Nyy(indexes);
                color_current_depth_range = obj_c(indexes ,:);
                active(t) = N;
                t= t+1;
                
                [y_run, x_run]= meshgrid((-(N-1)/2:(N-1)/2)*Hologram_sampling_interval,(-(N-1)/2:(N-1)/2)*Hologram_sampling_interval);
                r = sign(z)*sqrt(z^2 + y_run.^2 + x_run.^2);
                bb = exp(-1j*2*rand*pi)*exp(1j*k*r)./r;
                
                
                for index = 1: length(indexes)
                    
                    
                    pp = color_current_depth_range(index,color_index)*bb ;
                    
                    WRP(Nx(index):Nx(index)+N-1,Ny(index):Ny(index)+N-1) =  WRP(Nx(index):Nx(index)+N-1,Ny(index):Ny(index)+N-1)+ pp;
                    counter = counter + 1 ;
                    fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, total_wrp); %length(depth_ranges)
                    
                end
                
                
            end
            
            d_wrp_to_hologram = abs(max(d -  z_wrp));
            WRPHologram  = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
           
        end
        
        
        close all
        for p = 20:40 %reconstructed
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
        
        %
    end
    
    end_time = toc
    % imshow(rot90(uint8((original*2))),[]);
    %  uint8(round(RGB64*255));
    % figure;imshow((uint8((original))),[]);
    file_name = strcat('cup','_d_',num2str(d),'WRP_',num2str(length(depth_ranges)),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    % imwrite((uint8((original*2))),fullFileName);
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


% title('Segmentation in Conventional MWRP');
xlabel('No. of WRP')
ylabel('Calculation time(sec)')
hold on;
