function [ f ] = diffx2( psi )
    [nx,ny] = size(psi);
    hy = 2 * pi / ny;
    hx = hy;
    
    i=(1:nx)';
    j=(1:ny)';
    ip=circshift(i,-1);
    im=circshift(i, 1);
    
    f = (psi(ip,j)-psi(im,j)) / (2 * hx);    
end

