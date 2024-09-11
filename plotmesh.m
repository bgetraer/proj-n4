load('Models/Amundsen_Param.mat');
areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
el_length = sqrt(2*areas)*1E-3; % element side length in km
figure(1); clf; hold on;
patch( 'Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',el_length,'FaceColor','flat','EdgeColor','none')
axis tight equal
colormap(brewermap(100,'RdYlBu'))
set(gca,'yticklabel',get(gca,'ytick')*1E-3,'xticklabel',get(gca,'xtick')*1E-3,'fontsize',12)
xlabel('x (km)');
ylabel('y (km)');
cb = colorbar;
cb.Label.String = 'mean element side length (km)';
title('mesh resolution');

exportgraphics(gcf,'figures/mesh_map.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')

figure(2); clf; hold on;
hist(log10(el_length),100);
tick = [1:0.5:2.5, 5:5:25];
set(gca,'xtick',log10(tick));
set(gca,'xticklabel',tick,'fontsize',12)
xlabel('length (km)');
ylabel('count');
legend('mean element side length')
title('histogram of side lengths')
exportgraphics(gcf,'figures/mesh_hist.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
