function [div_tau] = smagorinsky_bilap_divergence(w,u,v,d,Cs)
    [alphax, alphay] = Smagorinsky_bilap_model(w, u, v, d);
    div_tau = Cs^4 * vector_divergence(alphax, alphay);
end

