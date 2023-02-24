function [ Cosp, k, dk ] = cospectrum_uv( sx, sy, lx, ly )
    [Cospx, ~,  ~] = cospectrum(sx,lx);
    [Cospy, k, dk] = cospectrum(sy,ly);
    Cosp = Cospx + Cospy;
end