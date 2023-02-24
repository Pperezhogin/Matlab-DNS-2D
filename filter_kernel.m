function [ uc ] = filter_kernel( u, hc, filter_type )
    % hc - filter width
    
    N = size(u,1);
    uf = fft2(u);
    kx = [0:N/2, (-N/2+1):-1];
    ky = [0:N/2, (-N/2+1):-1];
    [Kx, Ky] = ndgrid(kx, ky);    

    switch filter_type
        case 'hat_sum'
            G = (sinc(Kx * hc / 2 / pi) + sinc(Ky * hc / 2 / pi)) / 2;
        case 'hat_prod'
            G = sinc(Kx * hc / 2 / pi) .* sinc(Ky * hc / 2 / pi);    
        case 'spectral'
            G = abs((abs(Kx) < pi / hc).*(abs(Ky) < pi / hc));
        case 'gauss_sum'
            G = (exp(- Kx.^2 * hc^2 / 24) + exp(- Ky.^2 * hc^2 / 24)) / 2;
        case 'gauss_prod'
            G = exp(- (Kx.^2 + Ky.^2) * hc^2 / 24);            
    end
    
    uf = uf .* G;
    uc = ifft2(uf);    
end