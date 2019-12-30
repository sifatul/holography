%Main program of hologram generation by wavefront recording plane method
%by Openholo library project
%2017-10-30 update
%
addpath('D:\Mathlab\wrp\data');
addpath('D:\Mathlab\libs') ;
addpath('D:\Mathlab\wrp\data_scaled');
addpath('C:\Users\user\Desktop\WRP\trial1\good');
cd 'D:\My Research\holography\matlab\wrp'
clc;clear;close all;clear all; 


%% Input the object and prameter

three_object	
	
obj(:,1) = (obj(:,1)+0.0002).*2.5  ;	
obj(:,2) = obj(:,2)*2.2  ;	
obj(:,3) = obj(:,3)*80;	
 
obj_depth = max(obj(:,3)) - min(obj(:,3))
%%%%
% obj=openholo2(:,1:3);
lambda = [632.8e-9 532e-9 473e-9]; %RGB;                            % Wave length  
k = 2*pi./lambda;       
Hx=1280;                       % Hologram resolution  515   1025
Hy=720;
Hologram_sampling_interval = 8e-6;            % Hologram sampling interval 8e-6    3.9e-6


Cut=unique(obj(:,3));


%% setting WRP location
obj(:,3) = obj(1,3);
d = 0.15;
z=obj(:,3)+d;
size_obj = length(obj(:,1));
Hologram = zeros(Hy,Hx,3);
ROWS= Hy;                                     
COLS= Hx;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
[y_run, x_run]= meshgrid((-(Hx-1)/2:(Hx-1)/2)*Hologram_sampling_interval,(-(Hy-1)/2:(Hy-1)/2)*Hologram_sampling_interval);
tic
    for color_index = 1:length(lambda)
for o = 1: size_obj
    
    fprintf('%d\n',o);  
%     [y_run, x_run]= meshgrid((-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval,(-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval);
    yy=y_run-obj(o,1);
    xx=x_run-obj(o,2);
    r = sign(z(o))*sqrt(z(o)^2 + xx.^2 + yy.^2);
%     Sub_hologram = obj_c(o,1)*exp(1j*rand*2*pi)*exp(1j*k*r)./r;   

    Hologram(:,:,color_index)=Hologram(:,:,color_index)+obj_c(o,color_index)*exp(1j*rand*2*pi).*exp(sqrt(-1)*k(color_index)*r)./r;
%     figure,imshow(Hologram_wrp)
 
end
toc
% figure,imshow(angle(Hologram),[]);
%% Fresnel Propagation
end



%% reconstruction
for o=44  % reconstructed  
    for c=1:3
    close all;

        d2 = d+0.001 + o*0.001;
        original(:,:,c) = FresnelPropogation(k(c),v, h,-d2,Hologram(:,:,c));

    end
  figure; imshow(abs(original),[]);title(o);
end

phaseadd =WRPHologram;
phase_H = angle(phaseadd) + pi;
phase_H_image = uint8(255*phase_H/max(max(phase_H)));

FN=['wrp_hologram',num2str(Time),'_',num2str(count),'.bmp'];
imwrite(phase_H_image, FN, 'bmp');


