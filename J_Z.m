function [ f ] = J_Z( psi, w )
    % 
    f1 = diffx2(psi).*diffy2(w) - diffy2(psi).*diffx2(w); 
    f2 = diffy2(diffx2(psi).*w) - diffx2(diffy2(psi).*w);
    %f3 = diffx2(diffy2(w).*psi) - diffy2(diffx2(w).*psi);
    f = (f1 + f2) / 2;
end

