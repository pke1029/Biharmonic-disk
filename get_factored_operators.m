
function [r2, I_new, D1r_new, D2r_new, D3r_new, D4r_new] = get_factored_operators(rr, I, D1r, D2r, D3r, D4r)

r2 = (1-rr.^2);
I_new   = r2.*I;
D1r_new = r2.*D1r - 2*rr.*I;
D2r_new = r2.*D2r - 4*rr.*D1r - 2*I;
D3r_new = r2.*D3r - 6*rr.*D2r - 6*D1r;
D4r_new = r2.*D4r - 8*rr.*D3r - 12*D2r;

end