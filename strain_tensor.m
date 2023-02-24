function [ sxx, sxy, syy ] = strain_tensor( u, v )
    N = size(u,1);
    h = 2 * pi / N;
    
    idx = 1:N;
    idxp = circshift(idx,-1);
    idxm = circshift(idx,+1);    
    
    % in T points
    sxx = (u(idxp,:) - u) / h;
    syy = (v(:,idxp) - v) / h;
    % in w points
    sxy = (u - u(:,idxm)) / h + (v - v(idxm,:)) / h;
    
    sxy = sxy * 0.5;
end