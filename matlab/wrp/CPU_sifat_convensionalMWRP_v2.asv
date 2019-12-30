addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
cd 'D:\My Research\holography\matlab\wrp'
clc;clear;close all;clear all;



% % % zhao %maximum 120 wrp

%% medieval
% medieval
% obj(:,1) = (obj(:,1))*1.7;
% obj(:,2) = (obj(:,2))*1.7;
% obj(:,3) = obj(:,3)/7;


%%
% min_z = min(obj(:,2));
% max_z = max(obj(:,2));
% mid_z = (max_z-min_z)./2;
%
% obj1= obj(obj(:,2)<= min_z+mid_z,:);

%figure; plot3(obj1(:,1),obj1(:,2),obj1(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


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


% % % % four_items_new two_items
% % % % obj(:,1) = ((obj(:,1) -0.28)/300);
% % % % obj(:,2) = ((obj(:,2)+8  )/6200);
% % % % obj(:,3) = (obj(:,3)-0.5)/9;


%
% four_items_new
% obj(:,1) = ((obj(:,1)-0.15 )/800);
% obj(:,2) = ((obj(:,2)+2.5  )/1800);
% obj(:,3) = (obj(:,3)-0.5)/10;

% four_items
% four_items_new
% obj(:,1) = ((obj(:,1)-0.1 )/150);
% obj(:,2) = ((obj(:,2)+1.65  )/1800);
% obj(:,3) = (obj(:,3))/6.5;
% % zhao
% obj(:,1) = (obj(:,1)/750);
% obj(:,2) = (obj(:,2)/750);
% obj(:,3) = (obj(:,3)/6);
%


% zhao
% obj(:,1) = (obj(:,1)/350);
% obj(:,2) = (obj(:,2)/350);
% obj(:,3) = (obj(:,3)/4);


%%
% zhao
% obj(:,1) = (obj(:,1)/500);
% obj(:,2) = (obj(:,2)/500);
% obj(:,3) = (obj(:,3)/8);




% 1980*1024
% me_cube2
% obj(:,1) = ((obj(:,1)+0.2 )/800);
% obj(:,2) = ((obj(:,2)+0.15 )/600);
% obj(:,3) = (obj(:,3))/3;



% four_items_new
% obj(:,1) = ((obj(:,1)-0.1 )/130);
% obj(:,2) = ((obj(:,2)+1.65  )/1500);
% obj(:,3) = (obj(:,3))/6.5;
% papilon
% obj(:,1) = (obj(:,1)-0.0005)/2.5;
% obj(:,2) = (obj(:,2)-0.0006)/2.5;
% obj(:,3) = obj(:,3)/10 ;


% three_object
% %
% obj(:,1) = (obj(:,1)+0.0002).*4  ;
% obj(:,2) = (obj(:,2)-0.0002)*4  ;
% obj(:,3) = obj(:,3)*150;
% four_items_new
% obj(:,1) = ((obj(:,1)-0.16 )/300);
% obj(:,2) = (obj(:,2)+1.8 )/2500;
% obj(:,3) = (obj(:,3))/10;
% three_object
% 
% obj(:,1) = (obj(:,1)+0.0002).*2.5  ;
% obj(:,2) = obj(:,2)*2.2  ;
% obj(:,3) = obj(:,3)*80;


 



% temp  
% obj(:,1) = ((obj(:,1) )/300);
% obj(:,2) = (obj(:,2)+0.5 )/800;
% obj(:,3) = (obj(:,3))/12;


% three_object		
% 		
%  
% obj(:,1) = (obj(:,1)+0.0002).*3  ;
% obj(:,2) = obj(:,2)*3  ;
% obj(:,3) = obj(:,3)*150;
% % 		

% mm
% obj(:,1) = (obj(:,1)+1 )./780   ;
% obj(:,2) = obj(:,2)./620   ;
% obj(:,3) = obj(:,3)/6;
% d= 0.25; p = 166;

% books2		
% obj(:,1) = (obj(:,1)+0.5)./900   ;		
% obj(:,2) = (obj(:,2)+0.8)./2200   ;		
% obj(:,3) = obj(:,3)/8;		

% book_cube_chicken_cup			
% obj(:,1) = (obj(:,1)+0.5)./900   ;			
% obj(:,2) = (obj(:,2)-0.60)./350   ;			
% obj(:,3) = obj(:,3)/8;			

 
% me_cup
% obj(:,1) = (obj(:,1)+0.2)./600   ;
% obj(:,2) = (obj(:,2))./300   ;
% obj(:,3) = obj(:,3)/8; 
% 
	

		
		
% three_object		
% 		
% obj(:,1) = (obj(:,1)+0.0002).*2.5  ;		
% obj(:,2) = obj(:,2)*2.2  ;		
% obj(:,3) = obj(:,3)*80;		
		
		

% 		


% four_items  
% obj(:,1) = ((obj(:,1)-0.13 )/420);
% obj(:,2) = (obj(:,2) )/1200;
% obj(:,3) = (obj(:,3))/18;

% 
% three_object		
% obj(:,1) = (obj(:,1)+0.0002).*1.5  ;		
% obj(:,2) = obj(:,2)*1.5  ;		
% obj(:,3) = obj(:,3)*80;		
		
 

% two_cup
% obj(:,1) = (obj(:,1) )./350   ;
% obj(:,2) = (obj(:,2) )./450   ;
% obj(:,3) = obj(:,3)/8; 
% ironman		
% obj(:,1) = (obj(:,1)+0.25)./450   ;		
% obj(:,2) = (obj(:,2)-0.25)./350   ;		
% obj(:,3) = obj(:,3)/13;	 


% cube_cup_mug	
% obj(:,1) = (obj(:,1)+0.25)./800   ;	
% obj(:,2) = (obj(:,2)-0.15)./900   ;	
% obj(:,3) = obj(:,3)/12;

% ironman		
% obj(:,1) = (obj(:,1)+0.25)./400   ;		
% obj(:,2) = (obj(:,2)-0.40)./250   ;		
% obj(:,3) = obj(:,3)/13;	
% sad_life_increased
% obj(:,1) = (obj(:,1)-0.2)./350   ;     
% obj(:,2) = (obj(:,2)+2)./2500   ;     
% obj(:,3) = obj(:,3)/7 ; 
% d = 0.26;
% ironman_increased	
% obj(:,1) = (obj(:,1)+0.25)./600   ;     	
% obj(:,2) = (obj(:,2)-0.3)./500   ;     	
% obj(:,3) = obj(:,3)/7 ; 	

% papilon
% obj(:,1) = (obj(:,1))/3.5;
% obj(:,2) = (obj(:,2))/3.5;]
% obj(:,3) = obj(:,3)/8 ;
% three_object		
% 		
% obj(:,1) = (obj(:,1)+0.0002).*2.5  ;		
% obj(:,2) = obj(:,2)*2.2  ;		
% obj(:,3) = obj(:,3)*80;		

% % for 720X720
 
% 	
% % max(obj) - min(obj)	

% sad_life
% obj(:,1) = ((obj(:,1)-0.16 )/300);
% obj(:,2) = (obj(:,2)+2.25 )/2200;
% obj(:,3) = (obj(:,3))/20;
% 
% obj_depth = max(obj(:,3)) - min(obj(:,3))	
% d = 0.25;
	
four_items
obj(:,1) = ((obj(:,1)+0.15 )/350);
obj(:,2) = (obj(:,2) )/850;
obj(:,3) = (obj(:,3))/18;

 
 
max(obj)	
min(obj)	
% 	
% % max(obj) - min(obj)	
obj_depth = max(obj(:,3)) - min(obj(:,3))	
% d = 0.35;
d = 0.2;
obj_z = (obj(:,3));
[Cut,~,idx] = unique(obj_z);
obj_no = accumarray(idx(:),1);

figure;scatter([1:length(Cut)],Cut);
file_type = '.bmp';

%% Hologram parameter

Hologram_resolution_x = 1024;
Hologram_resolution_y = 1024;  % Hologram resolution
Hologram_resolution = strcat(num2str(Hologram_resolution_x),'X', num2str(Hologram_resolution_y)) ;                     % Hologram resolution
% Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval
Hologram_sampling_interval = 8e-6; %8e-6;%            % Hologram sampling interval

Nxx = ( round(obj(:,1)./Hologram_sampling_interval)+(Hologram_resolution_x)/2);
Nyy = (round (obj(:,2)./Hologram_sampling_interval)+(Hologram_resolution_y)/2 );



lambda = [632.8e-9 532e-9 473e-9]; %RGB;
k = 2*pi./lambda;

%% Fresnel propagation field

ROWS = Hologram_resolution_x;
COLS = Hologram_resolution_y;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));


%%







 
WRP_list= [2];
%   2, 4, 7, 14, 28, 49, 98
% WRP_list= [2, 4, 8,16,32,64];
prev_wrp_no = -1;
project_name = mfilename;
sub_dir = strcat('D:\Mathlab\wrp\final',obj_name,'',sprintf('_%0.9f',Hologram_sampling_interval),'\',num2str(Hologram_resolution),'\obj_depth',num2str(obj_depth),'\d_',num2str(d),'\layers_',sprintf('%d',length(Cut)),'\',project_name);
focus_location ='fly';
% mkdir(sub_dir);
% file_name = strcat('time', datestr(now, 'dd'),'.txt');
% time_file_name = fullfile(sub_dir, file_name);
% fileID = fopen(time_file_name,'a');

 
for ii= 1: length(WRP_list)
    
    original=zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    total_wrp = WRP_list(ii); 
    
    depth_ranges = depth_segmentation_conventional(Cut,total_wrp);
     
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
    
    layer_diff = abs(Cut(1)-Cut(2));
     
    WRPHologram = zeros(Hologram_resolution_x,Hologram_resolution_y,3);
    
    %
    tic
    for color_index = 1:length(lambda)
%         layer = 1;
       
        counter =1;
        for nwrp = 1:length(depth_ranges)
            WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);
            depth_range = cell2mat( depth_ranges(nwrp)) ;
            min_z = min(depth_range);
            max_z = max(depth_range);
            
            if mod(length(depth_range),2) == 1
                depth_range(end+1)= depth_range(end)+layer_diff;
                z_wrp = median(depth_range);
                depth_range(end) = [];
            else
                z_wrp = median(depth_range);  %  WRP location ;
            end
            indexes = find(obj_z>= min_z & obj_z <= max_z);
            
            z = (z_wrp-obj_z(indexes));
           
            
            
%                  N_all = round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1;  
%                     N_all = round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)));%sampling size of N
            %
            N_all = round(2*abs(z).*tan(lambda(color_index)./(2.*Hologram_sampling_interval))./Hologram_sampling_interval);
            
     
%             
%                         if min(N_all) ==0
%                             return
%                         end
            Nx = Nxx(indexes);
            Ny = Nyy(indexes);
            
          
            
            color_current_depth_range = obj_c(indexes ,color_index);
            
            
            for index = 1: length(indexes)
                
                [y_run, x_run]= meshgrid((-(N_all(index)-1)/2:(N_all(index)-1)/2)*Hologram_sampling_interval,(-(N_all(index)-1)/2:(N_all(index)-1)/2)*Hologram_sampling_interval);
                r = sign(z(index))*sqrt(z(index)^2 + y_run.^2 + x_run.^2);
                %                 pp = color_current_depth_range(index,color_index)*exp(-1j*2*rand*pi)*exp(1j*k*r)./r ;
                
                
                N=N_all(index);
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
                
                 Sub_hologram = color_current_depth_range(index)*exp(1j*rand*2*pi)*exp(1j*k(color_index)*r(1:xn,1:yn))./r(1:xn,1:yn) ;
                   

                WRP(Nx(index):nxn-1,Ny(index):nyn-1) =  WRP(Nx(index):nxn-1,Ny(index):nyn-1)+ Sub_hologram;
                %                 *exp(1j*rand*2*pi)
                %                 *exp(1j*k(c).*r(1:xn,1:yn))./r(1:xn,1:yn);
                
                %                 temp=zeros(Hologram_resolution_x+N(index), Hologram_resolution_y+N(index));
                %                 temp(Nx(index):Nx(index)+N(index)-1,Ny(index):Ny(index)+N(index)-1)=  Sub_hologram;
                %                 start_idx = round((N(index)+1)/2);
                %                 end_idx = round((N(index)-1)/2);
                %                 WRP=WRP+temp(start_idx:Hologram_resolution_x+end_idx,start_idx:Hologram_resolution_y+end_idx);
                
                %                 WRP(Nx(index):Nx(index)+N(index)-1,Ny(index):Ny(index)+N(index)-1) =  WRP(Nx(index):Nx(index)+N(index)-1,Ny(index):Ny(index)+N(index)-1)+ Sub_hologram;
                counter = counter + 1 ;
                fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, length(depth_ranges)); %length(depth_ranges)
                
            end
            
            d_wrp_to_hologram = abs(d -  z_wrp);
            WRPHologram(:,:,color_index)  = WRPHologram(:,:,color_index) + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
%             
        end
        
    end
    
     
    end_time = toc
      close all %120
    for p = 57
        for c= 1:3
            d2 = d+0.002 - p*0.0010
            original_red =  FresnelPropogation(k(c),v, h,-d2,gather(WRPHologram(:,:,c)));
            original_abs_red= abs((original_red));
            original(:,:,c)=255.*(original_abs_red./max(max(original_abs_red)));
            %             figure; imshow(rot90(abs(original(:,:,c)))*10,[]);title(p);
%                         toc;
%                                        file_name = strcat( 'color_',num2str(c),file_type);
%                         fullFileName = fullfile(sub_dir,file_name);
%                                 save_hologram(gather(WRPHologram(:,:,c)),fullFileName);
            %              imwrite((uint8((original))),fullFileName);
            %                       imwrite(uint8(rot90((original(:,:,co    `lor_index)))*10),fullFileName);
        end
                figure;   imshow(rot90(uint8((original*3))),[]);title(p);
    end
    

    %     fprintf('counter: %d color: %d wrp: %d \n',counter ,color_index, length(depth_ranges)); %length(depth_ranges)
    %     imshow(rot90(uint8((original*2))),[]);
    %  uint8(round(RGB64*255));
    %     figure;imshow((uint8((original))),[]);
    file_name = strcat('cube','_d_',num2str(d),'WRP_',num2str(length(depth_ranges)),'_total_layer_',num2str(length(Cut)),'_recon_d_',num2str(d2),file_type);
    fullFileName = fullfile(sub_dir,file_name);
    imwrite((uint8((original*2))),fullFileName);
    fprintf(fileID,' %d %0.3f \n ',length(depth_ranges), end_time);
    
end


%% reconstruction for red
fclose(fileID);
time = load(time_file_name);
x = time(:,1)
y = time(:,2)
sub_dir

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
