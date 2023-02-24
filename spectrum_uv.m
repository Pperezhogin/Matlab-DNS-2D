function [ Sp, k, dk ] = spectrum_uv( u,v )
    % spectrum on unit square, i.e. int(Sp) = int(u*u/2 + v*v/2) / L^2
    [Spu,~] = spectrum(u);
    [Spv,k,dk] = spectrum(v);
    Sp = Spu + Spv;
end