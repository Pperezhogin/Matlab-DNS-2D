function [w] = laplacian(psi)
    [u,v] = velocity(psi);
    w = vorticity(u,v);
end