function [ f ] = diffy( psi )
    [nx,ny] = size(psi);
    hy = 2 * pi / ny;
    
    i=(1:nx)';
    j=(1:ny)';
    jp=circshift(j,-1);
    jm=circshift(j, 1);
    
    f = (psi(i,jp) - psi(i,j)) / (hy);    
end

