function savePic(obj, Psi, s)
[~, z_max] = max(squeeze(sum(sum(abs(Psi).^2)) * obj.grid.dV));
t = s*obj.n_count_integrals*obj.dt;

cr = 10^6;

PsiIm = squeeze(Psi(:,:,z_max));
PhaseIm = squeeze(angle(Psi(:,:,z_max)));

f = figure('visible', 'off', 'Position', [10 10 900 400]);
p = uipanel('Parent',f,'BorderType','none'); 
p.Title = ['t = ', num2str(t, '%.2f'), ' s']; 
p.TitlePosition = 'centertop'; 
p.FontSize = 24;
          
subplot(1,2,1, 'Parent',p);
imagesc(cr*obj.grid.r, cr*obj.grid.r,abs(PsiIm).^2);
title(['\mid \Psi(x,y,', num2str(z_max, '%.0f'),') \mid^2']);

axis image;
ylabel('y, $\mu m$', 'interpreter', 'latex');
xlabel('x, $\mu m$', 'interpreter', 'latex');
set(gca,'YDir','normal')
colorbar;

subplot(1,2,2, 'Parent',p);
imagesc(cr*obj.grid.r, cr*obj.grid.r, PhaseIm);
title(['Arg(\Psi(x,y,', num2str(z_max, '%.0f'),'))']);
axis image;
ylabel('y, $\mu m$', 'interpreter', 'latex');
xlabel('x, $\mu m$', 'interpreter', 'latex');
set(gca,'YDir','normal')
colorbar;
caxis([-pi pi]);

saveas(f,[obj.Pictures,'Pic' num2str(s) '.png']);
clear f;
end