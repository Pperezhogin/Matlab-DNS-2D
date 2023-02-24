function [ w ] = vorticity( u, v )
    nx = size(u,1);
    ny = size(u,2);
    h = 2 * pi / ny;
    w = zeros(nx,ny);
    icx=[1:nx]';
    icmx=circshift(icx,1);
    icy=[1:ny]';
    icmy=circshift(icy,1);
    
    w(icx,icy)=(v(icx,icy)-v(icmx,icy)) / h - (u(icx,icy)-u(icx,icmy)) / h;
end
