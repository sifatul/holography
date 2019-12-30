
close all;
original_img = imread(' C:\Users\user\Desktop\journal\PSNR calculation\original_d_0.22WRP__total_layer_100_recon_d_0.166.bmp');
conv_img= imread('C:\Users\user\Desktop\journal\PSNR calculation\four_items_1024_1024\CPU_sifat_convensionalMWRP_v2\cube_d_0.2WRP_2_total_layer_100_recon_d_0.145.bmp');
proposed = imread('C:\Users\user\Desktop\journal\PSNR calculation\four_items_1024_1024\imid\cube_d_0.2WRP_3_total_layer_100_recon_d_0.14545.bmp');

% imshow(original_img);


crop =[540 470 100 100] %cube

original_cropped = imcrop(original_img,crop)*2;
%  crop =[260 400 100 200] %hen
I1 =imcrop(conv_img,crop);
% figure;
figure;imshow(I1);
figure; imshow(imcrop(proposed,crop));
figure; imshow(original_cropped)



peaksnr = psnr(imcrop(conv_img,crop),original_cropped)
peaksnr = psnr(imcrop(proposed,crop),original_cropped)
%
% addpath('D:\Mathlab\wrp\three_item_try_0.000008000\1980X1024\obj_depth0.083593\d_0.39\layers_96\imid');
% proposed= imread('book2EnhancedRMWRP_d_0.39threshold_0_total_layer_96_recon_d_0.315.bmp');
% figure;imshow(proposed);

I3 = imcrop(proposed,crop);

peaksnr = psnr(original_img,conv_img)
peaksnr = psnr(original_img,proposed)
%  figure;imshow(I3);

