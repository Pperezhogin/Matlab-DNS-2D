clc
clear all

nx = 128;
ny = 128;

hx = 2 * pi / nx;
hy = 2 * pi / ny;

xu = 0:hx:2*pi-hx;
yv = 0:hy:2*pi-hy;
yu = hy/2:hy:2*pi-hy/2;

[X,Y] = ndgrid(xu,yv);
[Xu,Yu] = ndgrid(xu,yu);

kd = 10;
U = 100;
dt = 0.5;
nu = 0.0001;

psi1 = randn(nx,ny) / 100000;
psi2 = randn(nx,ny) / 100000;
psi1 = psi1  - mean(psi1(:));
psi2 = psi2  - mean(psi2(:));

q1 = laplaceh(psi1) + 0.5*kd^2 * (psi2-psi1);
q2 = laplaceh(psi2) + 0.5*kd^2 * (psi1-psi2);

nt = 1;
t = 0;

while (1<2)
    q1p = q1;
    q2p = q2;
    q1 = q1 + dt/3 * (-J_Z(psi1,q1) + nu * laplaceh(q1) + U);
    q2 = q2 + dt/3 * (-J_Z(psi2,q2) + nu * laplaceh(q2) - U);
    [psi1,psi2] = find_psi(q1,q2,kd);
    
    q1 = q1p + dt/2 * (-J_Z(psi1,q1) + nu * laplaceh(q1) + U);
    q2 = q2p + dt/2 * (-J_Z(psi2,q2) + nu * laplaceh(q2) - U);
    [psi1,psi2] = find_psi(q1,q2,kd);
    
    q1 = q1p + dt * (-J_Z(psi1,q1) + nu * laplaceh(q1) + U);
    q2 = q2p + dt * (-J_Z(psi2,q2) + nu * laplaceh(q2) - U);
    [psi1,psi2] = find_psi(q1,q2,kd);
        
    pcolor(X,Y,laplaceh(psi1)); shading('flat')
    pause(0.01)
    
    t = t + dt;
    nt = nt + 1;
    t
end
