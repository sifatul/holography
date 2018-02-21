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
%% making quadrants quad1 and quad2

min_x = min(obj1(:,1));
max_x = max(obj1(:,1));
mid_x = (max_y-min_y)./2;

obj_quad1= obj1(obj1(:,1)<= min_x+mid_x,:); 
obj_quad2= obj1(obj1(:,1)> min_x+mid_x,:); 

%% setting WRP1 for obj_quad1

z_wrp = max(obj_quad1(:,2))+ 0.0002;  %  WRP location
z=z_wrp-obj_quad1(:,2);



N=round(abs(lambda.*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1;        %sampling size of N

Nx = round(obj_quad1(:,1)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;  
Ny = round(obj_quad1(:,3)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;

size_obj = length(obj_quad1(:,1));
Hologram_wrp1 = zeros(Hologram_resolution);

for o = 1: size_obj
     
    [y_run, x_run]= meshgrid((-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval,(-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval);
    r = sign(z(o))*sqrt(z(o)^2 + y_run.^2 + x_run.^2);
    Sub_hologram = exp(1j*rand*2*pi)*exp(1j*k*r)./r;   
    
    temp=zeros(Hologram_resolution+N(o), Hologram_resolution+N(o));    
    temp(Nx(o):Nx(o)+N(o)-1,Ny(o):Ny(o)+N(o)-1)= Sub_hologram;
    Hologram_wrp1=Hologram_wrp1+temp((N(o)+1)/2:Hologram_resolution+(N(o)-1)/2,(N(o)+1)/2:Hologram_resolution+(N(o)-1)/2);  
%     figure,imshow(Hologram_wrp)
 
end
%% reconstruct quad1
%{
ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
d=0.05;

%Hologram_wrp= Hologram_wrp2+ Hologram_wrp1;
WRPHologram1 = FresnelPropagation(Hologram_wrp1, Hologram_sampling_interval, Hologram_sampling_interval, d, lambda);

d2 = d+0.002 - 6*0.0002;
original = FresnelPropogation(k,v, h,-d2,WRPHologram1);
figure; imshow(abs(rot90(original)),[]);

%}
%% setting WRP2 for obj_quad2

z_wrp = max(obj_quad2(:,2))+ 0.0002;  %  WRP location
z=z_wrp-obj_quad2(:,2);


N=round(abs(lambda.*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1;        %sampling size of N

Nx = round(obj_quad2(:,1)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;  
Ny = round(obj_quad2(:,3)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;

size_obj = length(obj_quad2(:,1));
Hologram_wrp2 = zeros(Hologram_resolution);

for o = 1: size_obj
     
    [y_run, x_run]= meshgrid((-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval,(-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval);
    r = sign(z(o))*sqrt(z(o)^2 + y_run.^2 + x_run.^2);
    Sub_hologram = exp(1j*rand*2*pi)*exp(1j*k*r)./r;   
    
    temp=zeros(Hologram_resolution+N(o), Hologram_resolution+N(o));    
    temp(Nx(o):Nx(o)+N(o)-1,Ny(o):Ny(o)+N(o)-1)= Sub_hologram;
    Hologram_wrp2=Hologram_wrp2+temp((N(o)+1)/2:Hologram_resolution+(N(o)-1)/2,(N(o)+1)/2:Hologram_resolution+(N(o)-1)/2);  
%     figure,imshow(Hologram_wrp)
 
end

%% reconstruct quad2
%{
ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
d=0.05;

%Hologram_wrp= Hologram_wrp2+ Hologram_wrp1;
WRPHologram2 = FresnelPropagation(Hologram_wrp2, Hologram_sampling_interval, Hologram_sampling_interval, d, lambda);
d2 = d+0.002 - 6*0.0002;
original = FresnelPropogation(k,v, h,-d2,WRPHologram2);
figure;  imshow(abs(rot90(original)),[]);

%}

%% reconstruction in CGH

ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
d=0.05;

Hologram_wrp= Hologram_wrp2+ Hologram_wrp1;
WRPHologram = FresnelPropagation(Hologram_wrp, Hologram_sampling_interval, Hologram_sampling_interval, d, lambda);
d2 = d+0.002 - 6*0.0002;
original = FresnelPropogation(k,v, h,-d2,WRPHologram);
figure;  imshow(abs(rot90(original)),[]);


%{
    phaseadd =WRPHologram1;
    phase_H = angle(phaseadd) + pi;
    phase_H_image = uint8(255*phase_H/max(max(phase_H)));

    B = imrotate(phase_H_image,90);

    original_90 = FresnelPropogation(k,v, h,-d2,B);
    figure; imshow(rot90(abs(original_90)),[]);

    %imwrite(phase_H_image, 'wrp_hologram.bmp', 'bmp');
%}