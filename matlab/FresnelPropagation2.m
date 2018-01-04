function hologram = FresnelPropagation2(object, x, y, z, lambda)
   
    k = 2*pi/lambda;  
 
    H = exp(1i*k*z).*exp(-1i*pi*lambda*z*(x.^2+y.^2));   %Fourier transform of h
    O = fft2(object);                                   %Fourier transform of o
    hologram =ifft2(O.*H);          
 
end
