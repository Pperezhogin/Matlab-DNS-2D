function [ psi ] = inverse_laplaceh( w )
    N = size(w, 1);
    h = 2 * pi / N;
    kx = [0:N/2, (-N/2+1):-1];
    ky = [0:N/2, (-N/2+1):-1];
    [Kx, Ky] = ndgrid(kx, ky);
    Kx1 = sin(Kx*h/2)*2/h;
    Ky1 = sin(Ky*h/2)*2/h;
    iK2 = 1 ./ (Kx1.*Kx1 + Ky1.*Ky1);
    iK2(1,1) = 0;
    
    wf = fft2(w);
    psi = ifft2(-iK2 .* wf);    
end