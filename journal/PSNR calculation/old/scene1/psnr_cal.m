 
    close all;
    original_img = imread('C:\Users\user\Desktop\content\PSNR calculation\scene1\1280X720\original.bmp'); 
    proposed= imread('C:\Users\user\Desktop\content\PSNR calculation\scene1\1280X720\conv\hen_d_0.25WRP_8_total_layer_128_recon_d_0.173.bmp');
    conv_img = imread('C:\Users\user\Desktop\content\PSNR calculation\scene1\1280X720\imid\hen_d_0.25WRP_5_total_layer_128_recon_d_0.17329.bmp'); 

 imshow(conv_img);


 crop =[470 810 200 220] %cube 
%  crop =[260 400 100 200] %hen 
I1 =imcrop(conv_img,crop);
% figure;
figure;imshow(I1);
 

 
peaksnr = psnr(imcrop(conv_img,crop),imcrop(original_img,crop))
peaksnr = psnr(imcrop(proposed,crop),imcrop(original_img,crop))
% 
% addpath('D:\Mathlab\wrp\three_item_try_0.000008000\1980X1024\obj_depth0.083593\d_0.39\layers_96\imid');
% proposed= imread('book2EnhancedRMWRP_d_0.39threshold_0_total_layer_96_recon_d_0.315.bmp'); 
% figure;imshow(proposed);

I3 = imcrop(proposed,crop);

peaksnr = psnr(rmwrp_img,proposed)
%  figure;imshow(I3);

