
function [r, a, rr, aa, I, D1r, D2r, D3r, D4r, D2a, D4a, bc, bcinv] = get_operators(N, M)

if nargin == 0
    N = 25;
    M = 20;
end

% Chebyshev discretization in radial direction
N2 = (N-1)/2;
[D, r] = cheb(N);
r = r(r>0);
r = r';
I1 = eye(N2+1);

% Fourier discretization for azimuthal direction
M2 = M/2;
da = 2*pi/M; 
a = da*(1:M)';
D2a = toeplitz([-pi^2/(3*da^2)-1/6, .5*(-1).^(2:M)./sin(da*(1:M-1)/2).^2]);
D2a = kron(I1, D2a);
D4a = D2a*D2a;

% coordinates
[rr, aa] = meshgrid(r, a);
rr = rr(:);
aa = aa(:);

% derivative operators
D1r = D;
D2r = D1r*D;
D3r = D2r*D;
D4r = D3r*D;

% polar symmetry
I2 = eye(M);
I3 = [zeros(M2), eye(M2); eye(M2), zeros(M2)];
my_symm = @(D) kron(D(1:N2+1,1:N2+1), I2) + kron(D(1:N2+1,N+1:-1:N2+2), I3);
D1r = my_symm(D1r);
D2r = my_symm(D2r);
D3r = my_symm(D3r);
D4r = my_symm(D4r);

% identity operator
I = kron(I1, I2);

% boundary points
bc = @(D) D(M+1:end,M+1:end);
bcinv = @(u) [zeros(M,size(u,2)); u];

end