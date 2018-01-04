%% Test Point Cut ALL
close all;
clear all;
clc;

colorDevice = imaq.VideoDevice('kinect',1);
%Create a System object for the depth device.

depthDevice = imaq.VideoDevice('kinect',2);
%Initialize the camera.


%set up variables that don't need to be repeated
step(colorDevice);
step(depthDevice);

s = 1024; 
t = 1; 
lambda = 532e-9;
o=8;
k = 2*pi/lambda;
d = 0.25;
d2 = d - o*0.0001;
% Hologram_sampling_interval =6.4e-6;
Hologram_sampling_interval =7.4e-6;
                              
dx_cpu = Hologram_sampling_interval;   %      
dy_cpu = Hologram_sampling_interval;   %    
 
dx = gpuArray(dx_cpu);
dy = gpuArray(dy_cpu);


roi = [-1,0.5;-1,0.5;-1,1.3];

image_cpu = zeros(s);
film_cpu = zeros(s);
Hologram_cpu = zeros(s);


[Ny, Nx] = size(film_cpu); 
fx = 1./(Nx*dx_cpu);
fy = 1./(Ny*dy_cpu);  
x_cpu = ones(Ny,1)*[0:floor((Nx-1)/2) -ceil((Nx+1)/2)+1:-1]*fx;
y_cpu = [0:floor((Ny-1)/2) -ceil((Ny+1)/2)+1:-1]'*ones(1,Nx)*fy;

x=gpuArray(x_cpu);
y=gpuArray(y_cpu);
r=(x.^2+y.^2);
tic



%Acquire and view 5 frames of live Kinect point cloud data.

for i = 1:1   
    colorImage = step(colorDevice);  
    depthImage = step(depthDevice);

    ptCloud = pcfromkinect(depthDevice,depthImage,colorImage); 
    
    indices = findPointsInROI(ptCloud, roi);
    ptCloudB = select(ptCloud,indices);
    ptCloud=ptCloudB;
    
    [ptCloudOut,indices]= removeInvalidPoints(ptCloud); 

    Location = gpuArray(ptCloudOut.Location);


    % shift to positive axis 
    Obj(:,1) = Location(:,1) +min(Location(:,1))*-1 ; %x- axis
    Obj(:,2) = Location(:,2) +min(Location(:,2))*-1; %y- axis
    Obj(:,3) = Location(:,3) +min(Location(:,3))*-1; %z- axis




    Obj(:,4) = double (ptCloudOut.Color(:,1))./255 ; % R
    Obj(:,5) = double (ptCloudOut.Color(:,2))./255 ; % G
    Obj(:,6) = double (ptCloudOut.Color(:,3))./255 ; % B


    % scaling for 1024 
    Obj(:,1) = round (((s-10)/max (Obj(:,1))).* Obj(:,1) ) +1; % x- axis 
    Obj(:,2) =  round (((s-10)./max (Obj(:,2))).* Obj(:,2) ) +1; % y- axis 
    Obj(:,3) = round (((s-10)./max ( Obj(:,3) )).* Obj(:,3) ) +1; % z- axis 


    [C,ia,ic] = unique(Obj(:,1:2),'rows' , 'stable');
    Obj= Obj(ia,:);






    % Obj(:,1)= Obj(:,1)*t;
    % Obj(:,2)= Obj(:,2)*t;
    % Obj(:,3)= Obj(:,3)*1;


    % d = 0.5;
    % d  = Hologram_sampling_interval*(2*size_all)/lambda;
    %figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');ylabel ('z');title('Point Cloud');




    Obj = gather(Obj);


    z = Obj(:,3);
    A = categorical(z);
    Cut = str2double(categories(A));



    % % % GPU




    
    image = gpuArray(image_cpu);
    film = gpuArray(film_cpu);
    Hologram = gpuArray(Hologram_cpu);
    


    film1 = film;

    for i=1:length(Cut)
        O = Obj(z==Cut(i),:);
        O_image = film1; % for GPU
        
        test =sub2ind(size(O_image), O(:,1), O(:,2));
        O_image(test)=O(:,4);
        
        O_image = fft2(O_image); 
        d1 = d - Cut(i)*Hologram_sampling_interval/2;
        H = exp(1i*k*d1).*exp(-1i*pi*lambda*d1*r);   %Fourier transform of h
        
        film =O_image.*H;  
        film =ifft2(film);   
        
        Hologram = Hologram+film;
    end


    
    phase_H1 = angle(Hologram) + pi;
    phase_H_image = uint8(255*phase_H1/max(phase_H1)); %where is the use of phase_H_image

    dx = gather(dx);
    dy = gather(dy);
    Hologram = gather(Hologram);

toc

    originalR = FresnelPropagation2(Hologram, x,y, -d2, lambda);
    figure; imshow(abs((originalR)),[]);


end %end of realtime acquisition loop

release(colorDevice);
release(depthDevice);
