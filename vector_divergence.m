function [wim] = vector_divergence(lx ,ly)
    % lx, ly in v and u points
    N = size(lx,1);
    h = 2 * pi / N;
    
    idx = 1:N;
    idxm = circshift(idx,1);
    
    wim = - (lx - lx(idxm,:)) / h - (ly - ly(:,idxm)) / h;
end