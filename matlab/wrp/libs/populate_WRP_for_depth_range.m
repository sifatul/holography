function WRP = populate_WRP_for_depth_range(depth_range,Hologram_resolution_x,Hologram_resolution_y)
WRP = zeros(Hologram_resolution_x,Hologram_resolution_y);

z_wrp = median(depth_range);  %  WRP location ;
if  mod(length(depth_range),2) == 1
    z_wrp = z_wrp  +0.000005;
end
d_wrp_to_hologram = max(abs(d-z_wrp));
for layer = 1: length(depth_range)
    indexes = find(obj_z==depth_range(layer));
    z =   z_wrp - depth_range(layer);
    N =  (round(abs(lambda(color_index).*z./(Hologram_sampling_interval^2)/2)+0.5).*2-1 );        %sampling size of N
    Nx = Nxx(indexes);
    Ny = Nyy(indexes);
    color_current_depth_range = obj_c(indexes,color_index);
    [y_run, x_run]= meshgrid((-(N-1)/2:(N-1)/2)*Hologram_sampling_interval,(-(N-1)/2:(N-1)/2)*Hologram_sampling_interval);
    r =  ( sign(z)*sqrt(z^2 + y_run.^2 + x_run.^2));
    for n=1:length(indexes)
        counter = counter +1;
        fprintf('wrp %d counter: %d color: %d\n',length(depth_ranges), counter ,color_index);
        
        pp = color_current_depth_range(n)*exp(1j*rand*2*pi)*exp(1j*k*r)./r;
        
        
        WRP(Nx(n):Nx(n)+N-1,Ny(n):Ny(n)+N-1) =  WRP(Nx(n):Nx(n)+N-1,Ny(n):Ny(n)+N-1)+ pp ;
        
    end
    
end
         d_wrp_to_hologram = abs(max(d -  z_wrp));
            WRPHologram  = WRPHologram + FresnelPropagation((WRP), Hologram_sampling_interval, Hologram_sampling_interval, d_wrp_to_hologram, lambda(color_index));
            

end