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
lambdr = 532e-9;
lambdg = 532e-9;
lambdb = 473e-9;

o=8;

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
figure;

%above are constants 



%Acquire infinite frames
for P = 1:inf  

 tic
    colorImage = step(colorDevice);  
    depthImage = step(depthDevice);
 
    ptCloud = pcfromkinect(depthDevice,depthImage,colorImage); 
 
    indices = findPointsInROI(ptCloud, roi);
    ptCloudB = select(ptCloud,indices);
  
    ptCloud=ptCloudB;
     
    [ptCloudOut,indices]= removeInvalidPoints(ptCloudB); 
 
    Location = gpuArray(ptCloudOut.Location);
    
    Obj0=  gpuArray( Location );

    % shift to positive axis 
    Obj0(:,1) = Location(:,1) +min(Location(:,1))*-1 ; %x- axis
    Obj0(:,2) = Location(:,2) +min(Location(:,2))*-1; %y- axis
    Obj0(:,3) = Location(:,3) +min(Location(:,3))*-1; %z- axis



    Obj2=  gpuArray(ptCloudOut.Location );
    
    Obj2(:,1) = double (ptCloudOut.Color(:,1))./255 ; % R
    Obj2(:,2) = double (ptCloudOut.Color(:,2))./255 ; % G
    Obj2(:,3) = double (ptCloudOut.Color(:,3))./255 ; % B

    Obj=[gather(Obj0),gather(Obj2)];
    
    % scaling for 1024 
    Obj(:,1) = round (((s-10)/max (Obj(:,1))).* Obj(:,1) ) +1; % x- axis 
    Obj(:,2) =  round (((s-10)./max (Obj(:,2))).* Obj(:,2) ) +1; % y- axis 
    Obj(:,3) = round (((s-10)./max ( Obj(:,3) )).* Obj(:,3) ) +1; % z- axis 


    [C,ia,ic] = unique(Obj(:,1:2),'rows' , 'stable');
    Obj= Obj(ia,:);


    Obj = gather(Obj);


    z = Obj(:,3);
    A = categorical(z);
    Cut = str2double(categories(A));



    % % % GPU


    
    %image = gpuArray(image_cpu);
    film = gpuArray(film_cpu);
    Hologram_R = gpuArray(Hologram_cpu);
    Hologram_B = gpuArray(Hologram_cpu);
    Hologram_G = gpuArray(Hologram_cpu);


    film1 = film;

    for P=1:length(Cut)
        O = Obj(z==Cut(P),:);
        O_image_R = film1; % for GPU
        O_image_G = film1; % for GPU
        O_image_B = film1; % for GPU
        
        O_image_R(sub2ind(size(film1), O(:,1), O(:,2)))=O(:,4);
        O_image_G(sub2ind(size(film1), O(:,1), O(:,2)))=O(:,5);
        O_image_B(sub2ind(size(film1), O(:,1), O(:,2)))=O(:,6);
        
        O_image_R = fft2(O_image_R); 
        O_image_G = fft2(O_image_G); 
        O_image_B = fft2(O_image_B); 
        
        d1 = d - Cut(P)*Hologram_sampling_interval/2;
        H_R = exp(1i*(2*pi/lambdr)*d1).*exp(-1i*pi*lambdr*d1*r);   %Fourier transform of h
        H_G = exp(1i*(2*pi/lambdg)*d1).*exp(-1i*pi*lambdg*d1*r);   %Fourier transform of h
        H_B = exp(1i*(2*pi/lambdb)*d1).*exp(-1i*pi*lambdb*d1*r);   %Fourier transform of h
        
        film_R =O_image_R.*H_R;
        film_G =O_image_G.*H_G;
        film_B =O_image_B.*H_B;
        
        
        film_R =ifft2(film_R);
        film_G =ifft2(film_G);
        film_B =ifft2(film_B);
        
        Hologram_R = Hologram_R+film_R;
        Hologram_G = Hologram_G+film_G;
        Hologram_B = Hologram_B+film_B;
    end

    Hologram_R = gather(Hologram_R);
    Hologram_G = gather(Hologram_G);
    Hologram_B = gather(Hologram_B);

    originalR = FresnelPropagation2(Hologram_R, x,y, -d2, lambdr);
    originalG = FresnelPropagation2(Hologram_G, x,y, -d2, lambdg);
    originalB = FresnelPropagation2(Hologram_B, x,y, -d2, lambdb);
     %imshow(abs(rot90(originalR,-1)),[]);

      originalR= abs(rot90(originalR,-1));
      originalG= abs(rot90(originalG,-1));
      originalB= abs(rot90(originalB,-1));
     
      rgbImage = cat(3, originalR, originalG, originalB);
     
      imshow(rgbImage,[]);
      
      %RRR=abs(rot90(originalR,-1)); % storing red image 
     

toc
end %end of realtime acquisition loop

release(colorDevice);
release(depthDevice);
