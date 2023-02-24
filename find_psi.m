function [ psi1, psi2 ] = find_psi( q1, q2, kd )
    N = size(q1, 1);
    h = 2 * pi / N;
    kx = [0:N/2, (-N/2+1):-1];
    ky = [0:N/2, (-N/2+1):-1];
    [Kx, Ky] = ndgrid(kx, ky);
    Kx1 = sin(Kx*h/2)*2/h;
    Ky1 = sin(Ky*h/2)*2/h;
    K2 = Kx1.*Kx1 + Ky1.*Ky1;
    iK2 = 1 ./ (Kx1.*Kx1 + Ky1.*Ky1);
    iK2(1,1) = 0;
    iK2pkd2 = 1 ./ (Kx1.*Kx1 + Ky1.*Ky1+kd^2);
    
    q1f = fft2(q1);
    q2f = fft2(q2);
    
    idetA = iK2.*iK2pkd2;
    
    psi1f = - idetA .* (q1f.*(K2+0.5*kd^2) + q2f * 0.5*kd^2);
    psi2f = - idetA .* (q2f.*(K2+0.5*kd^2) + q1f * 0.5*kd^2);
    
    psi1 = ifft2(psi1f) - mean(q1(:)) / kd^2;
    psi2 = ifft2(psi2f) - mean(q2(:)) / kd^2;
    
    
end

