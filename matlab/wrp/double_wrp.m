%Main program of hologram generation by wavefront recording plane method
%by Openholo library project
%2017-10-30 update
%
clc;clear;close all;%% Input the object and prameter

load tea_pot.mat;
          % Hologram sampling interval

%figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

Obj(:,1)=(Obj(:,1)/40000 -0.0002 ).*1;
Obj(:,2)= (Obj(:,2)/20000 )+1.5;
Obj(:,3)= (Obj(:,3)/30000);

figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

%% 

lambda = 532e-9;                                % Wave length
k = 2*pi/lambda;       
Hologram_resolution=1025;                       % Hologram resolution     
Hologram_sampling_interval = 3.9e-6;            % Hologram sampling interval

min_z = min(Obj(:,2));
max_z = max(Obj(:,2));
mid_z = (max_z-min_z)./2;

obj1= Obj(Obj(:,2)<= min_z+mid_z,:); 

figure; plot3(obj1(:,1),obj1(:,2),obj1(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');

obj2= Obj(Obj(:,2) > min_z+mid_z,:); 
figure; plot3(obj2(:,1),obj2(:,2),obj2(:,3),'.');xlabel ('x');ylabel ('y');zlabel ('z');title('Point Cloud');


%% setting WRP1


z_wrp1 = max(obj1(:,2)) + 0.0002;  %  WRP location
z= z_wrp1-obj1(:,2);


N=round(abs(lambda.*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1;        %sampling size of N

Nx = round(obj1(:,1)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;  
Ny = round(obj1(:,3)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;

size_obj = length(obj1(:,1));
Hologram_wrp = zeros(Hologram_resolution);

for o = 1: size_obj
    
    %fprintf('%d\n',o);  
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

d= 0.05
WRPHologram1 = FresnelPropagation(Hologram_wrp, Hologram_sampling_interval, Hologram_sampling_interval, d, lambda);
d2 = d+0.002 - 6*0.0002;
original = FresnelPropogation(k,v, h,-d2,WRPHologram1);
figure; imshow(abs(original),[]);


%% setting WRP2

z_wrp2 = max(obj2(:,2)) + 0.0002;  %  WRP location
z = z_wrp2-obj2(:,2);

N=round(abs(lambda.*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1;        %sampling size of N

Nx = round(obj2(:,1)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;  
Ny = round(obj2(:,3)./Hologram_sampling_interval)+(Hologram_resolution-1)/2;

size_obj = length(obj2(:,1));
Hologram_wrp2 = zeros(Hologram_resolution);

for o = 1: size_obj
    
    %fprintf('%d\n',o);  
    [y_run, x_run]= meshgrid((-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval,(-(N(o)-1)/2:(N(o)-1)/2)*Hologram_sampling_interval);
    r = sign(z(o))*sqrt(z(o)^2 + y_run.^2 + x_run.^2);
    Sub_hologram = exp(1j*rand*2*pi)*exp(1j*k*r)./r;   
    
    temp=zeros(Hologram_resolution+N(o), Hologram_resolution+N(o));    
    temp(Nx(o):Nx(o)+N(o)-1,Ny(o):Ny(o)+N(o)-1)= Sub_hologram;
    Hologram_wrp2=Hologram_wrp2+temp((N(o)+1)/2:Hologram_resolution+(N(o)-1)/2,(N(o)+1)/2:Hologram_resolution+(N(o)-1)/2);  
%     figure,imshow(Hologram_wrp)
 
end

%% Fresnel Propagation

ROWS= Hologram_resolution;                                     
COLS= Hologram_resolution;
v=Hologram_sampling_interval.*(ones(COLS,1)*(-ROWS/2:ROWS/2-1))';
h=Hologram_sampling_interval.*(ones(ROWS,1)*(-COLS/2:COLS/2-1));
d= 0.05;

WRPHologram2 = FresnelPropagation(Hologram_wrp2, Hologram_sampling_interval, Hologram_sampling_interval, d, lambda);
  
d2 = d+0.002 - 6*0.0002;

original = FresnelPropogation(k,v, h,-d2,WRPHologram2);
figure; imshow(abs(original),[]);
%% final reconstruction

WRPHologram = WRPHologram2+ WRPHologram1;
d3 = 0.05; %0.081 any othe two CGS works best
original = FresnelPropogation(k,v, h,-d3,WRPHologram);
figure; imshow(rot90(abs(original)),[]);

%do interpolation or increase intensity in testing.m

