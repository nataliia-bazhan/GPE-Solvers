function showV(obj)
% show potential images from 2 sides

V = obj.getV();
          
% get dimensionless parameters
[l, E] = obj.getDimless();
r = obj.grid.r/l;
rz = obj.grid.rz/l;
Vxy = V(:, :, obj.grid.Nz/2)/E;
Vxz = squeeze(V(:, obj.grid.N/2, :))'/E;
          
figure('Position', [300 300 800 300])
          
subplot(1, 2, 1)
imagesc(r, r, Vxy)
axis square;
xlabel('x, $l_r$', 'interpreter', 'latex')
ylabel('y, $l_r$', 'interpreter', 'latex')
c = colorbar;
ylabel(c, 'V, $\hbar \omega_r$', 'interpreter', 'latex')
          
subplot(1, 2, 2)
imagesc(r, rz, Vxz)
axis square;
xlabel('x, $l_r$', 'interpreter', 'latex')
ylabel('z, $l_r$', 'interpreter', 'latex')
c = colorbar;
ylabel(c, 'V, $\hbar \omega_r$', 'interpreter', 'latex')
        
shg
end