function [ Z, k ] = spectrum_Z( w )
N = size(w,1);
L = 2*pi;
lmax = N/2;
kx = [0:N/2 (-N/2 + 1):-1];
ky = kx;
w1 = fft2(w);
Z = zeros(N/2-1,1);
for m=1:N
    for n=1:N
        K=sqrt( kx(n) * kx(n) + ky(m) * ky(m) );
        l=floor(K);
        if (l>0)&&(l<lmax)
            Z(l)=Z(l) + w1(n,m)*conj(w1(n,m));
        end
    end
end
Z = Z * L^2 / (2 * N^4);
k = 1.5:(N/2-0.5);
k=k';
end