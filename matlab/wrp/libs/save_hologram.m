function save_hologram( hologram,file_name )
 

phase_H1 = angle(hologram) + pi;
phase_H_image = uint8(255*phase_H1/max(max(phase_H1)));
% figure; imshow(abs(phase_H_image),[]);
imwrite(phase_H_image,file_name); 
end

