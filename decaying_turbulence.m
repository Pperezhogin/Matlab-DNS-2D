clc
clear all
restoredefaultpath;

nx = 128;
ny = 128;

hx = 2 * pi / nx;
hy = 2 * pi / ny;

xu = 0:hx:2*pi-hx;
yv = 0:hy:2*pi-hy;
yu = hy/2:hy:2*pi-hy/2;

[X,Y] = ndgrid(xu,yv);
[Xu,Yu] = ndgrid(xu,yu);

psi = randn(nx,ny);
psi = psi - mean(psi(:));
psi = filter_kernel(psi,2*pi/10,'gauss_prod');
w = laplacian(psi);
[u,v] = velocity(psi);
maxu = max(abs(u(:)));

nt = 1;
t = 0;
dt = hx / maxu / 3;
Cs = 0.1;
bf_width = hx * sqrt(6);

while (t<33)
    rw1 = compute_rhs(w,bf_width,Cs);
    w1 = w + rw1 * dt / 2;
    rw2 = compute_rhs(w1,bf_width,Cs);
    w2 = w + rw2 * dt/2;
    rw3 = compute_rhs(w2,bf_width,Cs);
    w3 = w + rw3 * dt;
    rw4 = compute_rhs(w3,bf_width,Cs);
    
    phi = dt/6 * (rw1 + 2 * rw2 + 2 * rw3 + rw4);
    
    inv_phi = inverse_laplaceh(phi);
    %gamma = - 2 * dot(w(:),inv_phi(:)) / dot(inv_phi(:), phi(:));
    gamma = 1;
    
    w = w + gamma * phi;    
    
    if (mod(nt,10) == 0)
        %imagesc(w);
        %pause(0.01)
        [u,v] = velocity(inverse_laplaceh(w));
        current_CFL = gamma * max(max(abs(u(:))), max(abs(v(:)))) * dt / hx
        gamma
        E = mean(mean(u.^2+v.^2))/2
        imagesc(w)
        pause(0.05)
        t
    end 
        
    t = t + dt;
    nt = nt + 1;
end