function [ Sp, k, dk ] = spectrum( w )
    % spectrum on unit square, i.e. int(Sp) = int(w*w/2) / L^2
    [Sp,k,dk] = cospectrum(w,w/2);
    Sp = max(Sp,0);
    N = size(w,1);
    lmax = max(find(k<N/2+0.5));
    Sp = Sp(1:lmax);    
    k = k(1:lmax);
    dk = dk(1:lmax);
end