function [ Cosp, k, dk ] = cospectrum( w, tau )
    %cospectrum on unit square, i.e. int(Cosp) = int(w*tau) / L^2
%     N = size(w,1);
%     lmax = floor(sqrt((N/2)^2 + (N/2)^2));
%     kx = [0:N/2 (-N/2 + 1):-1];
%     ky = kx;
%     wf = fft2(w);
%     tauf = fft2(tau);
%     Cosp = zeros(lmax,1);
%     
%     for m=1:N
%         for n=1:N
%             K=sqrt( kx(n) * kx(n) + ky(m) * ky(m) );
%             l=floor(K);
%             if (l>0)
%                 Cosp(l)=Cosp(l) + real(conj(wf(n,m)) * tauf(n,m));
%             end
%         end
%     end
% 
%     Cosp = Cosp / (N^4);
%     k = 1:1:lmax;
%     k=k';
%     dk = ones(size(k));

    N = size(w,1);
    lmax = floor(sqrt((N/2)^2 + (N/2)^2));
    kx = [0:N/2 (-N/2 + 1):-1];
    ky = kx;
    wf = fft2(w);
    tauf = fft2(tau);
    
    k_cusp = 10; % wavenumber from which reduce density of points
    k = 1:k_cusp; k = k';
    dk = ones(k_cusp-1,1);
    for i = k_cusp+1:lmax
        k_new = k_cusp/(k_cusp-1)*k(i-1);
        if ( k_new > lmax )
            dk(i-1) = dk(i-2);
            break;
        else
           k(i) = k_new; 
           dk(i-1) = k(i) - k(i-1);
        end
    end
    
   lmax = length(k); 
   Cosp = zeros(lmax,1);
   
   Cosp_matrix = real(conj(wf).*tauf);
   Cosp_matrix(1,1) = 0; % to omit constant
   
   [Kx,Ky] = ndgrid(kx,ky);
   K_matrix = sqrt(Kx.^2+Ky.^2);
   numerator = log((K_matrix+1e-15) / k_cusp);
   denominator = log(k_cusp / (k_cusp-1));
   mask = K_matrix < k_cusp - 1e-15;
   idx = floor(K_matrix).*mask + (k_cusp + floor(numerator/denominator)).*(1-mask);
   idx(1,1) = 1;
      
    for m=1:N
        for n=1:N
            Cosp(idx(n,m))=Cosp(idx(n,m)) + Cosp_matrix(n,m);
        end
    end

    Cosp = Cosp / (N^4) ./ dk; 
end