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
                              
dx = Hologram_sampling_interval;   %      
dy = Hologram_sampling_interval;   %    
 


roi = [-1,0.5;-1,0.5;-1,1.3];

image = zeros(s);
film = zeros(s);
Hologram = zeros(s);


[Ny, Nx] = size(film); 
fx = 1./(Nx*dx);
fy = 1./(Ny*dy);  
x = ones(Ny,1)*[0:floor((Nx-1)/2) -ceil((Nx+1)/2)+1:-1]*fx;
y = [0:floor((Ny-1)/2) -ceil((Ny+1)/2)+1:-1]'*ones(1,Nx)*fy;

r=(x.^2+y.^2);
figure;

%above are constants 



%Acquire infinite frames
for i = 1:inf  

    colorImage = step(colorDevice);  
    depthImage = step(depthDevice);
 
    ptCloud = pcfromkinect(depthDevice,depthImage,colorImage); 
 
    indices = findPointsInROI(ptCloud, roi);
    ptCloudB = select(ptCloud,indices);
  
    ptCloud=ptCloudB;
     
    [ptCloudOut,indices]= removeInvalidPoints(ptCloudB); 
 
    Location = (ptCloudOut.Location);
 
    
    %Obj0= Location
    % shift to positive axis 
    Obj(:,1) = double (Location(:,1) +min(Location(:,1))*-1) ; %x- axis
    Obj(:,2) = double (Location(:,2) +min(Location(:,2))*-1); %y- axis
    Obj(:,3) = double (Location(:,3) +min(Location(:,3))*-1); %z- axis



   %Obj2=  gpuArray(ptCloudOut.Location );
    
    Obj(:,4) = double (ptCloudOut.Color(:,1))./255 ; % R
    Obj(:,5) = double (ptCloudOut.Color(:,2))./255 ; % G
    Obj(:,6) = double (ptCloudOut.Color(:,3))./255 ; % B


    
    % scaling for 1024 
    Obj(:,1) = round (((s-10)/max (Obj(:,1))).* Obj(:,1) ) +1; % x- axis 
    Obj(:,2) =  round (((s-10)./max (Obj(:,2))).* Obj(:,2) ) +1; % y- axis 
    Obj(:,3) = round(((s-10)./max ( Obj(:,3) )).* Obj(:,3) ) +1; % z- axis 


    [C,ia,ic] = unique(Obj(:,1:2),'rows' , 'stable');
    Obj= Obj(ia,:);


    figure; plot3(Obj(:,1),Obj(:,2),Obj(:,3),'.');xlabel ('x');ylabel ('y');ylabel ('z');title('Point Cloud');




    %Obj = gather(Obj);


    z = Obj(:,3);
    A = categorical(z);
    Cut = str2double (categories(A));


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


    
    phase_H = angle(Hologram) + pi;
    phase_H_image = uint8(255*phase_H/max(max(phase_H)));
   % imwrite(phase_H_image, 'testing.bmp', 'bmp');
   
    clear Obj;
    clear Hologram;
    clear phase_H_image;
   % Hologram = gather(Hologram);



    originalR = FresnelPropagation2(Hologram, x,y, -d2, lambda);
    imshow(abs(rot90(originalR,-1)),[]);


end %end of realtime acquisition loop

release(colorDevice);
release(depthDevice);
