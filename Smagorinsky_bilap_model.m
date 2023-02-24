function [alphax, alphay] = Smagorinsky_bilap_model(w, u, v, d)
    N = size(w,1);
    h = 2 * pi / N;
    
    idx = 1:N;
    idxp = circshift(idx,-1);
    idxm = circshift(idx,+1);
    
    [ sxx, sxy, syy ] = strain_tensor(u, v);
    s2w = sxy.^2;
    
    nuT = d^4*sqrt(2*(sxx.^2 + syy.^2 + 2*(s2w+s2w(idxp,:)+s2w(idxp,idxp)+s2w(:,idxp))/4));
    nux = (nuT + nuT(:,idxm)) / 2;
    nuy = (nuT + nuT(idxm,:)) / 2;
    
    wx = (w(idxp,:) - w) / h;
    wy = (w(:,idxp) - w) / h;
    
    
    alphax = nux.*laplacian(wx);
    alphay = nuy.*laplacian(wy);
end