function [wn, psin, un, vn] = update_state(w, psi, u, v, rhs_w, rhs_psi, rhs_u, rhs_v, dt)
    wn = w + dt * rhs_w;
    psin = psi + dt * rhs_psi;
    un = u + dt * rhs_u;
    vn = v + dt * rhs_v;
end