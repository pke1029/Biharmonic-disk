
function my_plot(r, a, zz)

N = length(r);
M = length(a);

% add boundary points
[aa, rr] = meshgrid([0;a], [r,-r(end)]);
Z = reshape(zz, [M,N]);
Z = [Z, circshift(Z(:,end), M/2)];
Z = [Z(end,:); Z];
% mesh([r,-r(end)], [0;a], Z)

% interpolate to fine grid
aq = linspace(0, 2*pi, 150)';
rq = 0:0.04:1;
[aaq, rrq] = meshgrid(aq, rq);
Zq = interp2(aa, rr, Z', aaq, rrq, 'spline');
% mesh(aaq, rrq, Zq)

% convert to cartesian coordinates
[xx, yy] = pol2cart(aaq, rrq);
surf(xx, yy, real(Zq), 'EdgeColor', 'none', 'FaceColor', 'interp')
axis square

end