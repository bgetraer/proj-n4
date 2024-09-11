
loaddata=0;
if loaddata % {{{
	disp('loading domain data')

	% BedMachine data	
	xLim = [-2E6, -1E6];
	yLim = [-1E6,0.2E6];
	x=linspace(xLim(1),xLim(2),5E3);
	y=linspace(yLim(1),yLim(2),5E3);
	[X Y]=meshgrid(x,y);
	bed = interpBedmachineAntarctica(X,Y,'bed');
	mask = interpBedmachineAntarctica(X,Y,'mask');
	rock = mask==1;

	% outline of ice domain in x (m)
	xExp=[-1669315.4719493026,-1669315.4719493026,-1193987.0047960179,...
		-1026979.7055259449,-1026979.7055259449,-1556906.7128252152,...
		-1772089.1945770399,-1772089.1945770399,-1669315.4719493026];
	% outline of ice domain in y (m)
	yExp=[-420940.0927398402,-829553.2314259715,-829553.2314259715,...
		-530867.1000391102,-58750.3117179424,170008.8123696489,...
		70446.7685740285,-420940.0927398402,-420940.0927398402];

	ADD_coastline=shaperead('../../add_coastline_medium_res_line_v7_9/add_coastline_medium_res_line_v7_9.shp');
	ADD_rockpolygon=shaperead('../../add_rock_outcrop_medium_res_polygon_v7/add_rock_outcrop_medium_res_polygon_v7.3.shp');
	IBCSO_shelfbreakline=shaperead('../../Shelf_break_ANTARCTICA_DAmblas/Shelf_break_Antarctica_DAmblas.shp');

	x1=[];x2=[];y1=[];y2=[];
	for i=1:length(ADD_coastline)
		switch ADD_coastline(i).surface
			case {'ice rumples','grounding line'}
				x1=[x1, ADD_coastline(i).X];
				y1=[y1, ADD_coastline(i).Y];
			otherwise 
				x2=[x2, ADD_coastline(i).X];
				y2=[y2, ADD_coastline(i).Y];
		end
	end
	xrock=[];yrock=[];
	for i=1:length(ADD_rockpolygon)
		xrock=[xrock, ADD_rockpolygon(i).X];
		yrock=[yrock, ADD_rockpolygon(i).Y];
	end
	xshelf=IBCSO_shelfbreakline.X(1:end-1);
	yshelf=IBCSO_shelfbreakline.Y(1:end-1);
	save('domainplot','xExp','yExp','x1','y1','x2','y2','xrock','yrock','xshelf','yshelf');
else
	load('domainplot')
end
% }}}


figure(1);clf;
set(gcf,'position',get(gcf,'position').*[1 1 0 0] + [0 0 1120 850])

ax1=axes; hold on;
h_bed=imagesc(x,y,bed);
set(ax1,'ydir','normal')
axis equal off

ax2=axes; hold on;
set(ax2,'color','none');
shade=hillshade(imgaussfilt(bed,10),x,y,'azimuth',150,'altitude',30)./256;
h_shade=imagesc(x,y,shade,'alphadata',shade.^3-0.45);
h_rockimage = imagesc(x,y,~rock,'alphadata',rock);
h_exp=plot(xExp,yExp,'-','color',[1,0.4,0.4],'linewidth',3);
h_gl=plot(x1,y1,'-k','linewidth',0.75);
h_coast=plot(x2,y2,'-k','linewidth',2);
plot(xrock,yrock,'-k','linewidth',1);
h_rock=plot([nan nan nan],[nan nan nan],'sk','markerfacecolor','k');
%h_shelf=plot(xshelf,yshelf,'--k','linewidth',2);
h_shelf=plot([nan nan nan],[nan nan nan],'--k','linewidth',0.75);

axis equal
set([ax1,ax2],'xlim',xLim,'ylim',yLim)
set(ax2,'yticklabel',get(gca,'ytick')*1E-3,'xticklabel',get(gca,'xtick')*1E-3,'fontsize',18)
xlabel('x (km)');
ylabel('y (km)');


%yline(0);xline(0);
hlgd=legend(ax2,[h_exp h_coast h_gl h_rock h_shelf], 'model domain','coastline/ice front','grounding line',...
	'rock outcrop','continental shelf break','Location','northeastoutside');
ax1.Position=ax2.Position;


cmap1=flip(brewermap(1000,'Spectral'));
cmap2=gray;
set(ax1,'colormap',cmap1,'clim',[-2500,3750]);
set(ax2,'colormap',cmap2,'clim',[0 1]);

cb=colorbar(ax1,'horiz','axislocation','in','Ticks',[-2500:1250:3750],'Position',[hlgd.Position(1)+0.01,hlgd.Position(2)-0.1,hlgd.Position(3)-0.02,0.04]);
set(cb,'FontSize',10,'TickLength',0.02);
title(cb,'Bed topography (m)','FontSize',14);


ax3=axes;
hold on;
h_shelf=patch(xshelf,yshelf,[1,1,1],'linestyle','--');
h_exp=patch(xExp,yExp,[1,0.4,0.4],'edgecolor','k');
h_gl=plot(x1,y1,'-k','linewidth',0.75);
h_coast=plot(x2,y2,'-k','linewidth',2);
%plot(xrock,yrock,'.k','markersize',0.001);
h_rock=plot([nan nan nan],[nan nan nan],'sk','markerfacecolor','k');

axis equal tight
set(gca,'xtick',[],'ytick',[],'box','off','linewidth',1)
set(gca,'yticklabel',{},'xticklabel',{},'color','none','YColor','none','XColor','none')
set(gca,'fontsize',18)

ax3.Position=[hlgd.Position(1)-0.055,0.11,0.4,0.4];

%% printing
printfig=0;
if printfig
exportgraphics(gcf,'domain.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
set([ax1.Children;ax2.Children;ax3.Children;hlgd;cb;ax2],'Visible','off')
set([h_bed,h_shade,h_rockimage],'Visible','on')
exportgraphics(gcf,'domainimage.png','BackgroundColor','none','Resolution',600,'ContentType','image')
end
