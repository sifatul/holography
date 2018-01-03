%% Test Point Cut ALL
close all;
clear all;
clc;



colorDevice = imaq.VideoDevice('kinect',1);
%Create a System object for the depth device.

depthDevice = imaq.VideoDevice('kinect',2);
%Initialize the camera.

step(colorDevice);
step(depthDevice);
%Load one frame from the device.

colorImage = step(colorDevice);
depthImage = step(depthDevice);
%Extract the point cloud.

ptCloud = pcfromkinect(depthDevice,depthImage,colorImage);
%Initialize a point cloud player to visualize 3-D point cloud data. The axis is set appropriately to visualize the point cloud from Kinect.





%Acquire and view 5 frames of live Kinect point cloud data.

for i = 1:5    
   colorImage = step(colorDevice);  
   depthImage = step(depthDevice);
 
   ptCloud = pcfromkinect(depthDevice,depthImage,colorImage); 
   [ptCloudOut,indices]= removeInvalidPoints(ptCloud); 
  
     %colormap(gray)
   % pcshow(ptCloud)
   %view(player,ptCloudOut);
end


min_of_obj = abs((min(min(ptCloudOut.Location)))) ;
Obj(:,1) = ptCloudOut.Location(:,1)+min_of_obj ; %x- axis
Obj(:,2) = ptCloudOut.Location(:,2)+min_of_obj ; %y- axis
Obj(:,3) = ptCloudOut.Location(:,3) +min_of_obj ; %z- axis




Obj(:,4) = double (ptCloudOut.Color(:,1))./255 ; % R
Obj(:,5) = double (ptCloudOut.Color(:,2))./255 ; % G
Obj(:,6) = double (ptCloudOut.Color(:,3))./255 ; % B

%pcwrite(ptCloud,'sameple_pointCloud2','PLYFormat','binary');
%fileID = fopen('test_pcl.txt','w');


%fprintf(fileID,' %8f %8f %8f %8f %8f %8f \n', [ptCloudOut.Location,double(ptCloudOut.Color)]);

%fclose(fileID);



% Obj = load('cube15c.txt');
%Obj = load('Tang25.txt');
% Obj = load('mi25.txt');
% Obj = load('3obj.txt');


%{
Obj(:,1)= round((Obj(:,1)+100).*5); %15
Obj(:,2)= round((Obj(:,2)+50).*6); %15
Obj(:,3) = -(Obj(:,3)-230);  % tang25
%}

%[~,ind,~]=unique(Obj(:, 1:2:3), 'rows', 'stable');
%Obj=Obj(sort(ind), :);



s = 1024; % size
t = 1; % 放大倍数，以512为基准


Obj(:,1) = round ((s./max (Obj(:,1))).* Obj(:,1) ) -2; % x- axis 
Obj(:,2) =  round ((s./max (Obj(:,2))).* Obj(:,2) -2); % y- axis 
Obj(:,3) = round ((s./max ( Obj(:,3) )).* Obj(:,3) -2); % z- axis 


[C,ia,ic] = unique(Obj(:,1:2),'rows' , 'stable');
Obj= Obj(ia,:);


lambda = 532e-9;
o=8;

k = 2*pi/lambda;
% Hologram_sampling_interval =6.4e-6;
Hologram_sampling_interval =7.4e-6;
                              
dx = Hologram_sampling_interval;   %      
dy = Hologram_sampling_interval;   %    

Obj(:,1)= Obj(:,1)*t;
Obj(:,2)= Obj(:,2)*t;
Obj(:,3)= Obj(:,3)*1;

d = 0.25;
% d = 0.5;
% d  = Hologram_sampling_interval*(2*size_all)/lambda;
figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');ylabel ('z');title('Point Cloud');


release(colorDevice);
release(depthDevice);



image = zeros(s);
film = zeros(s);
Hologram = zeros(s);
z = Obj(:,3);
A = categorical(z);
Cut = str2double(categories(A));

[Ny, Nx] = size(film); 
fx = 1./(Nx*dx);
fy = 1./(Ny*dy);  
x = ones(Ny,1)*[0:floor((Nx-1)/2) -ceil((Nx+1)/2)+1:-1]*fx;
y = [0:floor((Ny-1)/2) -ceil((Ny+1)/2)+1:-1]'*ones(1,Nx)*fy;

% % % GPU



dx = gpuArray(dx);
dy = gpuArray(dy);
image = gpuArray(image);
film = gpuArray(film);
Hologram = gpuArray(Hologram);
x=gpuArray(x);
y=gpuArray(y);


film1 = film;



tic
for i=1:length(Cut)
    O = Obj(z==Cut(i),:);
    O_image = film1; % for GPU
%     O_image = zeros(s); % for CPU
    test =sub2ind(size(O_image), O(:,1), O(:,2));
    O_image(test)=O(:,4);
%      O_image(sub2ind(size(O_image), O(:,1), O(:,2)))=O(:,4);
%     image = image + O_image;
        O_image = fft2(O_image); 
    d1 = d - Cut(i)*Hologram_sampling_interval/2;
    H = exp(1i*k*d1).*exp(-1i*pi*lambda*d1*(x.^2+y.^2));   %Fourier transform of h
           
    %Fourier transform of o
     film =O_image.*H;  
     
    film =ifft2(film);          

%     film = FresnelPropagation(O_image, dx, dy, d1, lambda);
    Hologram = Hologram+film;
end
toc
dx = gather(dx);
dy = gather(dy);
Hologram = gather(Hologram);


phase_H1 = angle(Hologram) + pi;
phase_H_image = uint8(255*phase_H1/max(phase_H1)); %where is the use of phase_H_image


    d2 = d - o*0.0001;
% %     original = zeros(s);
    originalR = FresnelPropagation2(Hologram, x,y, -d2, lambda);
    figure; imshow(abs(rot90(originalR)),[]);


