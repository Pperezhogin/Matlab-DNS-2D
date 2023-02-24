function [u,v] = velocity(psi)
    N = size(psi,1);
    h = 2 * pi / N;
    
    idx = 1:N;
    idxp = circshift(idx,-1);
    
    u = - (psi(:,idxp) - psi) / h;
    v = + (psi(idxp,:) - psi) / h;
end