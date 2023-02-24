function [rhs_w] = compute_rhs(w,hx,Cs)
    psi = inverse_laplaceh(w);
    [u,v] = velocity(psi);
    rhs_w = -J_Z(psi,w) + smagorinsky_bilap_divergence(w,u,v,hx,Cs);
end

