function hologram = FresnelPropagation2(object, x, y, z, lambda)
   
    k = 2*pi/lambda;
    %[Ny, Nx] = size(object);

    %temp_x = [1:Nx]-floor((Nx+1)/2);
    %temp_y = [1:Ny]-floor((Ny+1)/2);
    
    
    %fx = 1./(Nx*dx);
    %fy = 1./(Ny*dy);  
    %x = ones(Nx,1)*[0:floor((Nx-1)/2) -ceil((Nx+1)/2)+1:-1]*fx;          %Note order of points for FFT
   % y = [0:floor((Ny-1)/2) -ceil((Ny+1)/2)+1:-1]'*ones(1,Nx)*fy;
 
    H = exp(1i*k*z).*exp(-1i*pi*lambda*z*(x.^2+y.^2));   %Fourier transform of h
    O = fft2(object);                                   %Fourier transform of o
    hologram =ifft2(O.*H);          
   % max_phase_step = max(pi*lambda*z*(x(1, 1:(Nx-1)).^2 - x(1, 2:(Nx)).^2));
    %du = dx;
    %dv = dy;
end
