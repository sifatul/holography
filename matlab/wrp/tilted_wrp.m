clc;clear;close all;%% Input the object and prameter

load tea_pot.mat;
          % Hologram sampling interval

%figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

Obj(:,1)=(Obj(:,1)/40000 -0.0002 ).*1;
Obj(:,2)= (Obj(:,2)/20000 )+1.5;
Obj(:,3)= (Obj(:,3)/30000);
lambda = 532e-9;                                % Wave length
k = 2*pi/lambda;       
Hologram_resolution=1025;                       % Hologram resolution     
Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval




%figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');
%% making object 1
min_y = min(Obj(:,2));
max_y = max(Obj(:,2));
mid_y = (max_y-min_y)./2;

obj1= Obj(Obj(:,2)<= min_y+mid_y,:); 

figure;plot3(obj1(:,1),obj1(:,2),obj1(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

%% setting WRP1

wrp_x = obj1(:,1);
wrp_z = obj1(:,3);
%wrp_y = obj1(:,2).*0 + max(obj1(:,2))+0.0002; theta = 0; file_name='tilted_0.jpg';
%wrp_y = obj1(:,2).*0 + max(obj1(:,2))+ 0.0238; theta = pi/18; file_name='tilted_10.jpg';
wrp_y = obj1(:,2).*0 + max(obj1(:,2))+0.6238; theta = pi/4; file_name='tilted_45.jpg';
%wrp_y = obj1(:,2).*0 + max(obj1(:,2))+1.4988; theta = pi/3; file_name='tilted_60.jpg';
%figure;plot3(wrp_x, wrp_y,wrp_z ,'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');



X = wrp_x;
Y = wrp_y*cos(theta) - wrp_z*sin(theta);
Z = wrp_y*sin(theta) + wrp_z*cos(theta);

%figure;plot3(X, Y,Z ,'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

z= Y-obj1(:,2);


N=round(abs(lambda.*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1;        %sampling size of N

Nx = round(obj1(:,1)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;  
Ny = round(obj1(:,3)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;

size_obj = length(obj1(:,1));
Hologram_wrp = zeros(Hologram_resolution);

for o = 1: size_obj
     
    [y_run, x_run]= meshgrid((-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval,(-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval);
    r = sign(z(o))*sqrt(z(o)^2 + y_run.^2 + x_run.^2);
    Sub_hologram = exp(1j*rand*2*pi)*exp(1j*k*r)./r;   
    
    temp=zeros(Hologram_resolution+N(o), Hologram_resolution+N(o));    
    temp(Nx(o):Nx(o)+N(o)-1,Ny(o):Ny(o)+N(o)-1)= Sub_hologram;
    Hologram_wrp=Hologram_wrp+temp((N(o)+1)/2:Hologram_resolution+(N(o)-1)/2,(N(o)+1)/2:Hologram_resolution+(N(o)-1)/2);  
%     figure,imshow(Hologram_wrp)
 
end

%% Fresnel Propagation

ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
d=0.05;
%d= z_wrp1*300;
WRPHologram1 = FresnelPropagation(Hologram_wrp, Hologram_sampling_interval, Hologram_sampling_interval, d, lambda);





for o=0:1:20  % reconstructed  
    d2 = d+0.002 - o*0.0002;
    original = FresnelPropogation(k,v, h,-d2,WRPHologram1);
    figure; imshow(abs(original),[]);
end
d2 = d+0.002 - 9*0.0002; %0 degree
original = FresnelPropogation(k,v, h,-d2,WRPHologram1);
h= figure; imshow(abs(rot90(original)),[]);

saveas(h,file_name);


%{
phaseadd =WRPHologram1;
phase_H = angle(phaseadd) + pi;
phase_H_image = uint8(255*phase_H/max(max(phase_H)));

B = imrotate(phase_H_image,90);

original_90 = FresnelPropogation(k,v, h,-d2,B);
figure; imshow(rot90(abs(original_90)),[]);

%imwrite(phase_H_image, 'wrp_hologram.bmp', 'bmp');
%}