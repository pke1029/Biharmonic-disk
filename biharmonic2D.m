
clear

N = 35;
M = 20;
[r, a, rr, aa, I, D1r, D2r, D3r, D4r, D2a, D4a] = get_operators(N, M);
[r2, I, D1r, D2r, D3r, D4r] = get_factored_operators(rr, I, D1r, D2r, D3r, D4r);

% biharmonic operator
L = D4r + 2*D2r*D2a./rr.^2 + D4a./rr.^4 + 2*D3r./rr - 2*D1r*D2a./rr.^3 - D2r./rr.^2 + 4*D2a./rr.^4 + D1r./rr.^3;
L = L * diag(1./r2);
% boundary condition
L = L(M+1:end,M+1:end);

% compute eigenvalues
[V, lam] = eig(L); 
lam = diag(lam); 
[lam,ii] = sort(lam);
V = V(:,ii);
V = [zeros(M,size(V,2)); V];
V = V.*r2;

% plot
for i = 1:9
    subplot(3,3,i)
    my_plot(r, a, V(:,i))
    axis off
    view([0,90])
    title(lam(i))
end