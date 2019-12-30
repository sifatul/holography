addpath('D:\Mathlab\wrp\finalthree_obj_0.000008000\1980X1024\obj_depth0.07228\psnr');


rmwrp = rgb2gray(imread('rmwrp.bmp'));
conv= rgb2gray(imread('conventional.bmp'));
imid =rgb2gray( imread('imid.bmp'));

corr2(conv,rmwrp)
corr2(imid,rmwrp)


addpath('C:\Users\user\Desktop\papilion');

rmwrp = rgb2gray(imread('rmwrp.bmp'));
conv= rgb2gray(imread('conven.bmp'));
imid =rgb2gray( imread('imid.bmp'));



cor = corr2(rmwrp,conv)
cor = corr2(rmwrp,imid)

% peaksnr = psnr(A,ref)

% rmwrp_img = (imread('rmwrp.bmp'));
% conv_img= (imread('conven.bmp'));
% imid_img =( imread('imid.bmp'));
% 
% 
% 
% % cor = corr2(rmwrp,conv)
% % cor = corr2(rmwrp,imid)
% 
% crop =[260 250 200 200]; %% [xmin ymin width height]
% rmwrp = imcrop(rmwrp_img,crop);
% 
%  
% conv = imcrop(conv_img,crop);
% imid = imcrop(imid_img,crop);
% 
% % figure;
% imshow(rmwrp);
% 
% peaksnr = psnr(conv,rmwrp)
% peaksnr = psnr(imid,rmwrp)