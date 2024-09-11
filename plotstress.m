clc;
disp('Running plotting routine');
readdata=1;
if readdata
clear variables;
	fname='Models/Amundsen_n3_ccsm85_TransientRun.mat';
	md=loadmodel(fname);
	y2s=365*24*60*60; % s/yr
	vx=cell2mat({md.results.TransientSolution.Vx}); % array of vx on vertices (m/y)
	vy=cell2mat({md.results.TransientSolution.Vy}); % array of vy on vertices (m/y)

	t=cell2mat({md.results.TransientSolution.time}); % time (y)
end


% prepare data to analyze
vx0 = md.initialization.vx; % initial velocity (x component) (m/y)
vy0 = md.initialization.vy; % initial velocity (y component) (m/y)
t1=2013+15; % second time to compare to
vx1=vx(:,find(t==t1)); % velocity at t1 (x component) (m/y)
vy1=vy(:,find(t==t1)); % velocity at t1 (y component) (m/y)
% ocean mask
mask=(md.initialization.vx==0 & md.initialization.vy==0);
vx0(mask) = NaN;          % do not consider non-moving domain 
vy0(mask) = NaN;          % do not consider non-moving domain
vx1(mask) = NaN;     % do not consider non-moving domain
vy1(mask) = NaN;     % do not consider non-moving domain

% get effective strain rates from initialization velocity
sr0 = strainrate_SSA(md,vx0,vy0, 0); % strain rates per element (1/y)
e0  = sr0.eff;								 % effective strain rates per element (1/y)
le0 = log10(e0);                       % log of effective strain rates 
caxmin = round(nanmean(le0) - 2*nanstd(le0)); 
caxmax = ceil(nanmean(le0) + 3*nanstd(le0));

% total mass lost
vaf0=md.results.TransientSolution(1).IceVolumeAboveFloatation; % m^3
vaf1=md.results.TransientSolution(find(t==t1)).IceVolumeAboveFloatation; % m^3
dmaf=md.materials.rho_ice.*(vaf1-vaf0); % kg


% standard colormap
cax=[-3,0];
cm=brewermap(17,'spectral');
cmpos=flip(cm(1:round(end/2),:));
cmneg=cm(round(end/2):end,:);
%% FIGURE 1: DISTRIBUTION OF STRAIN RATE VALUES {{{
%figure(1);clf;hold on;
%hist(le0,100);
%N=hist(le0,100);
%xline([caxmin caxmax],'linewidth',1.5);
%% mark off std and mean
%for i=0:3
%	if i==0
%		xline(nanmean(le0),'-r','linewidth',1.5);
%		text(nanmean(le0),max(N)*3/4,'\mu','color','r')
%	else
%		xline(nanmean(le0)-i*nanstd(le0),'-r','linewidth',1.5);
%		xline(nanmean(le0)+i*nanstd(le0),'-r','linewidth',1.5);
%		text(nanmean(le0)-i*nanstd(le0),max(N)*3/4,sprintf('-%i\\sigma',i),'color','r')
%		text(nanmean(le0)+i*nanstd(le0),max(N)*3/4,sprintf('+%i\\sigma',i),'color','r')
%	end
%end
%text(caxmin-0.2,max(N)/2,'min caxis threshold','Rotation',90,'HorizontalAlignment','left')
%text(caxmax+0.2,max(N)/2,'max caxis threshold','Rotation',90,'HorizontalAlignment','left')
%xlabel('log(effstrainrate)')
%ylabel('number of elements')
%% }}}
% FIGURE 2: MAP OF INITIAL EFFECTIVE STRAIN RATES {{{
fig2=figure(2);clf;
ax1=axes;
patch( 'Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',le0,'FaceColor','flat','EdgeColor','none')
axis equal tight
% colorbar
set(ax1,'colormap',cmpos,'clim',cax);
cb=colorbar('ticks',cax(1):1/3:cax(2));
ticklab=cell(1,length(cb.Ticks));
for i=1:length(cb.Ticks)
	if mod(cb.Ticks(i),1)==0
		ticklab{i}=['10^{',num2str(cb.Ticks(i)),'}'];
	end
end
set(cb,'TickLabels',ticklab)
title(cb,'$(a^{-1})$','Interpreter','latex','FontSize',14);

% geometric axes
set(ax1,'xticklabels',get(ax1,'xtick').*1E-3,'yticklabels',get(ax1,'ytick').*1E-3,'fontsize',12);
xlabel('UPS grid easting (km)'); ylabel('UPS grid northing (km)');
title('${\dot{\varepsilon}}_{e}$','interpreter','latex','fontsize',16)

set(ax1, 'Layer', 'top','Color','none')
%exportgraphics(fig2,'e_eff.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
% }}}

% plot the change in effective strain rate
sr1 = strainrate_SSA(md,vx1,vy1, 0); % strain rates per element (1/y)
e1 = sr1.eff;            % effective strain rates per element (1/y)
d = (e1-e0); % difference in effective strain rates per element (1/y)
d(d==0)=NaN; 
ld = log10(abs(d)); % log difference in effective strain rates per element (1/y)

%% FIGURE 3: DISTRIBUTION OF CHANGE IN EFFECTIVE STRAIN RATE VALUES {{{
%% positive and negative values of change in eff. strain rate IN LOG SPACE
%caxmin = -4;
%caxmax = 0;
%yLIM=ylim;
%figure(3);clf;hold on;
%hneg=histogram(ld(d<0),100);
%hpos=histogram(ld(d>0),100);
%xline(nanmean(ld),'r','linewidth',1.5)
%text(nanmean(ld),yLIM(2)*31/32,'\mu','color','r')
%xline([caxmin,caxmax],'r','linewidth',1.5)
%text(caxmin-0.2,yLIM(2)/2,'min caxis threshold','Rotation',90,'HorizontalAlignment','left')
%text(caxmax+0.2,yLIM(2)/2,'max caxis threshold','Rotation',90,'HorizontalAlignment','left')
%lg=legend([hneg hpos],'$\Delta\dot\varepsilon < 0 $','$\Delta\dot\varepsilon > 0 $',...
%	'interpreter','latex', 'location','nw');
%xlabel('log$_{10}(|\Delta\dot\varepsilon|)$','interpreter','latex')
%ylabel('number of elements')
%% }}}
%% FIGURE 4: MAP OF CHANGE IN EFFECTIVE STRAIN RATES {{{
%figure(4);clf;hold on;
%plotmodel(md,'data',ld,'caxis',[caxmin caxmax],'figure',4);
%% colorbar ticks and label
%cb=get(gca,'colorbar');
%ticklab=cell(1,length(cb.Ticks));
%for i=1:length(cb.Ticks)
%	if cb.Ticks(i)>-3
%		ticklab{i}=num2str(10.^cb.Ticks(i));
%	else
%		ticklab{i}=['10^{',num2str(cb.Ticks(i)),'}'];
%	end
%end
%set(cb,'TickLabels',ticklab)
%cb.Label.Interpreter = 'latex';
%cb.Label.String = '$|\Delta\dot\varepsilon|$ $(a^{-1})$';
%% axes
%set(gca,'xticklabels',get(gca,'xtick').*1E-3,'yticklabels',get(gca,'ytick').*1E-3);
%xlabel('UPS grid easting (km)');
%ylabel('UPS grid northing (km)');
%colormap(brewermap(15,'spectral'))
%% }}}
% FIGURE 4: MAP OF CHANGE IN EFFECTIVE STRAIN RATES {{{
fig4=figure(4);clf;
ld_pos=nan(size(ld)); ld_neg=nan(size(ld));
ld_pos(d>0)=ld(d>0);  ld_neg(d<0)=ld(d<0);
ld_pos((e1-e0)==0)=-100; ld_neg((e1-e0)==0)=-100;
ax1=axes;
patch( 'Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',ld_neg,'FaceColor','flat','EdgeColor','none')
axis equal tight
ax2=axes;
patch( 'Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',ld_pos,'FaceColor','flat','EdgeColor','none')
axis equal tight
linkaxes([ax1,ax2])
%%Hide the top axes
ax1.Visible = 'off'; ax1.XTick = []; ax1.YTick = [];
% COLORMAPS
colormap(ax1,cmneg)
colormap(ax2,cmpos)

set([ax1,ax2],'Position',[.05 .11 .685 .815],'CLim',cax,'fontsize',12);
%cb1 = colorbar(ax1,'Position',[.75 .11 .0675 .815]);
%cb2 = colorbar(ax2,'Position',[.825 .11 .0675 .815]);
cb1 = colorbar(ax1,'Position',[0.6375    0.1098                  0.0381    0.8154/2-0.005],'ydir','reverse');
cb2 = colorbar(ax2,'Position',[0.6375    0.1098+0.8154/2+0.005   0.0381    0.8154/2-0.005]);

% colorbar ticks and label
ticklabpos=cell(1,length(cb1.Ticks));
ticklabneg=cell(1,length(cb1.Ticks));
for i=1:length(cb1.Ticks)
	if mod(cb1.Ticks(i),1)==0
		ticklabpos{i}=['+10^{',num2str(cb1.Ticks(i)),'}'];
		ticklabneg{i}=['-10^{',num2str(cb1.Ticks(i)),'}'];
	end
end
ticklabneg{1}='';ticklabpos{1}=['\pm10^{',num2str(cb1.Ticks(1)),'}'];
title(cb2,'$(a^{-1})$','Interpreter','latex','Fontsize',14);
set(cb1,'TickLabels',ticklabneg);
set(cb2,'TickLabels',ticklabpos);
% axes
set(ax2,'xticklabels',get(ax2,'xtick').*1E-3,'yticklabels',get(ax2,'ytick').*1E-3);
xlabel(ax2,'UPS grid easting (km)'); ylabel(ax2,'UPS grid northing (km)');
title(ax2,'$\Delta\dot\varepsilon_e$','Interpreter','latex','fontsize',16)
set(ax2, 'Layer', 'top','Color','none')
%exportgraphics(fig4,'delta_e_eff.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
% }}}


adresults=load('./Models/gradient.mat','gradient'); % gradient of VAF wrt n
g=md.materials.rho_ice.*adresults.gradient; % gradient of MAF wrt n
g(g==0)=NaN;
ng=g;
lg = log10(abs(ng));
% FIGURE 5: MAP OF  GRAD dM/dn {{{
fig5=figure(5);clf;
lg_pos=nan(size(lg)); lg_neg=nan(size(lg));
lg_pos(g>0)=lg(g>0);  lg_neg(g<0)=lg(g<0);
lg_pos((e1-e0)==0)=-100; lg_neg((e1-e0)==0)=-100;
ax1=axes;
patch( 'Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',lg_neg,'FaceColor','flat','EdgeColor','none')
axis equal tight
ax2=axes;
patch( 'Faces', md.mesh.elements, 'Vertices', [md.mesh.x md.mesh.y],'CData',lg_pos,'FaceColor','flat','EdgeColor','none')
axis equal tight
linkaxes([ax1,ax2])
%%Hide the top axes
set(ax1,'Visible','off','XTick',[],'YTick', []); set(ax2, 'Layer', 'top','Color','none') 
% COLORMAPS
caxg=[4 9]+3;
cm=brewermap(17,'spectral');
cmpos=flip(cm(1:round(end/2),:));
cmneg=cm(round(end/2):end,:);
colormap(ax1,cmneg)
colormap(ax2,cmpos)
set([ax1,ax2],'Position',[.05 .11 .685 .815],'CLim',caxg,'fontsize',12);

cb1 = colorbar(ax1,'Position',[0.6375    0.1098                  0.0381    0.8154/2-0.005],'ydir','reverse');
cb2 = colorbar(ax2,'Position',[0.6375    0.1098+0.8154/2+0.005   0.0381    0.8154/2-0.005]);
% colorbar ticks and label
ticklabpos=cell(1,length(cb1.Ticks));
ticklabneg=cell(1,length(cb1.Ticks));
for i=1:length(cb1.Ticks)
   if mod(cb1.Ticks(i),1)==0
      ticklabpos{i}=['+10^{',num2str(cb1.Ticks(i)),'}'];
      ticklabneg{i}=['-10^{',num2str(cb1.Ticks(i)),'}'];
   end
end
ticklabneg{1}='';ticklabpos{1}=['\pm10^{',num2str(cb1.Ticks(1)),'}'];
title(cb2,'$kg/m^2$','Interpreter','latex','Fontsize',14);
set(cb1,'TickLabels',ticklabneg);
set(cb2,'TickLabels',ticklabpos);
% axes
set(ax2,'xticklabels',get(ax2,'xtick').*1E-3,'yticklabels',get(ax2,'ytick').*1E-3);
xlabel(ax2,'UPS grid easting (km)'); ylabel(ax2,'UPS grid northing (km)');
%title(ax2,'$\frac{1}{|\Delta\mathrm{MAF}|}\frac{\partial \mathrm{MAF}}{\partial n}$','Interpreter','latex','fontsize',16)
title(ax2,'$\mathcal{D}M(n)$','Interpreter','latex','fontsize',16)
%exportgraphics(fig5,'dMdn.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
%exportgraphics(fig5,'dMdn_normalized.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
% }}}

% normalize and compare the three
en=normalize(le0);
dn=normalize(ld);
gn=normalize(lg);

axmin=floor(min([en;dn;gn]));
axmax=ceil(max([en;dn;gn]));

figure(6);clf;
subplot(1,3,1);hold on;
xlim([axmin axmax])
ylim([axmin axmax])
line([xlim],[ylim],'color','k')
plot(en,dn,'.');
axis square
r2 = 1 - nanvar(dn-en)/nanvar(dn);
title(sprintf('R^2 = %0.2f',r2))

subplot(1,3,2);hold on;
xlim([axmin axmax])
ylim([axmin axmax])
line([xlim],[ylim],'color','k')
plot(en(g<0),gn(g<0),'.');
plot(en(g>0),gn(g>0),'.');
axis square
r2 = 1 - nanvar(gn-en)/nanvar(gn);
title(sprintf('R^2 = %0.2f',r2))

subplot(1,3,3);hold on;
xlim([axmin axmax])
ylim([axmin axmax])
line([xlim],[ylim],'color','k')
plot(dn(g<0),gn(g<0),'.');
plot(dn(g>0),gn(g>0),'.');
axis square
r2 = 1 - nanvar(gn-dn)/nanvar(gn);
title(sprintf('R^2 = %0.2f',r2))

error();


%%%%
figcrossplot=figure(10);clf; hold on
x=d(~isnan(d)&~isnan(g));
y=g(~isnan(d)&~isnan(g));

C=-10
lx = sign(x).*(log10(1+abs(x)/(10^C)));
ly = sign(y).*(log10(1+abs(y)/(10^C)));
histogram2(lx,ly,200,'DisplayStyle','tile');
m = mean(abs(ly./lx));
plot(xlim,m*xlim,'--k');
plot(xlim,-m*xlim,'--k');
xtick=sort([10.^[-9:3:0] -10.^[-9:3:0]]);
lxtick=sign(xtick).*(log10(1+abs(xtick)/(10^C)));
set(gca,'xtick',sort([lxtick 0]),'xticklabels',sort([xtick 0]))

%%%%
figcrossplot=figure(11);clf; hold on
x=e0(~isnan(d)&~isnan(g));
y=g(~isnan(d)&~isnan(g));

C=-10
lx = sign(x).*(log10(1+abs(x)/(10^C)));
ly = sign(y).*(log10(1+abs(y)/(10^C)));
histogram2(lx,ly,200,'DisplayStyle','tile');
m = mean(abs(ly./lx));
plot(xlim,m*xlim,'--k');
plot(xlim,-m*xlim,'--k');




% FIGURE 6: COMPARISON OF EFF STR RATE TO CHANGE IN EFF STRAIN RATE {{{
figure(6);clf;hold on;
plot(e0,d,'.')
xx=[0 10.^(-10:0.2:0)];
threshold=0.35;
hmax=plot(xx,threshold*ones(size(xx)),'-r','linewidth',1); % 
hmin=plot(xx,-threshold*ones(size(xx)),'-r','linewidth',1); % 
h1000=plot(xx,xx*1E1,'-k'); % 1000% increase: most fit tightly within this
h100=plot(xx,xx*1E0,'-k'); % 100% increase
hp1=plot(xx,xx*1E-3,'-k'); % 0.1% increase
hnp1=plot(xx,-xx*1E-3,'-k'); % 0.1% decrease
hn100=plot(xx,-xx*1E0,'-k'); % 100% decrease: mathematical limit);
symlog(gca, 'x', -5);
symlog(gca, 'y', -10);

text(0.5,hmax.YData(end),'0.35 a^{-1}','color','r','verticalalignment','bottom') % 1% increase
text(0.5,hmin.YData(end),'-0.35 a^{-1}','color','r','verticalalignment','top') % 1% increase

text(h1000.XData(end),h1000.YData(end),'1000% increase') % 1% increase
text(h100.XData(end),h100.YData(end),'100% increase') % 1% increase
text(hp1.XData(end),hp1.YData(end),'0.1% increase') % 1% increase
text(hnp1.XData(end),hnp1.YData(end),'0.1% decrease') % 1% increase
text(hn100.XData(end),hn100.YData(end),'100% decrease') % 1% increase

xlabel('$\dot\varepsilon$ $(a^{-1})$','interpreter','latex')
ylabel('$\Delta \dot\varepsilon$ $(a^{-1})$','interpreter','latex')
ytick=get(gca,'YTick');
yticklab=get(gca,'YTickLabel');
set(gca,'YTick',ytick(1:2:end),'YTickLabel',{yticklab{1:2:end}})
% }}}
% FIGURE 7: COMPARISON OF EFF STR RATE TO GRADIENT {{{
figure(7);clf;hold on;
plot(e0,g,'.')
symlog(gca, 'x', -6);
symlog(gca, 'y', -6);
% }}}
% FIGURE 8: COMPARISON OF EFF STR RATE TO GRADIENT {{{
figure(8);clf;
subplot(1,6,2:5);hold on;
hall=plot(d,g,'.')
hthresh=plot(d(abs(g)>1E7),g(abs(g)>1E7),'.');
symlog(gca, 'x', -9);
symlog(gca, 'y', -8);
x=hall.XData;
y=hall.YData;
ytick=get(gca,'YTick');
yticklab=get(gca,'YTickLabel');
set(gca,'YTick',ytick(1:2:end),'YTickLabel',{yticklab{1:2:end}},'fontsize',8)
yl=ylim;
[N_neg,b_neg] = hist(y(x<0),100);
[N_pos,b_pos] = hist(y(x>0),100);
yln=min(abs(hthresh.YData));
yline([-yln,yln])
ytick=get(gca,'ytick');
yticklab=get(gca,'yticklabel');
grid minor; grid on;
xlabel('$\Delta \dot\varepsilon$ $(a^{-1})$','interpreter','latex')

subplot(1,6,1);hold on;
barh(b_neg,N_neg);
barh(b_neg((abs(b_neg)>yln)),N_neg(abs(b_neg)>yln))
yline([-yln,yln])
xlim([0,max([N_neg,N_pos])])
ylim(yl)
set(gca,'XDir','reverse','ytick',ytick,'yticklabel',yticklab)
grid on
ylabel('gradient of VAF wrt n')
subplot(1,6,6);hold on;
barh(b_pos,N_pos);
barh(b_pos((abs(b_pos)>yln)),N_pos(abs(b_pos)>yln))
yline([-yln,yln])
xlim([0,max([N_neg,N_pos])])
ylim(yl)
set(gca,'ytick',ytick,'yticklabel',yticklab)
grid on
% }}}



x=hall.XData;
y=hall.YData;
histogram(y(x>0),100)
histogram(y(x<0),100)

% FIGURE 8: COMPARISON OF EFF STR RATE TO GRADIENT {{{
figure(8);clf;hold on;
plot(d(abs(g)>1E6),g(abs(g)>1E6),'.')
symlog(gca, 'x', -9);
symlog(gca, 'y', 5);
% }}}



text(0.5,hmax.YData(end),'0.35 a^{-1}','color','r','verticalalignment','bottom') % 1% increase
text(0.5,hmin.YData(end),'-0.35 a^{-1}','color','r','verticalalignment','top') % 1% increase

text(h1000.XData(end),h1000.YData(end),'1000% increase') % 1% increase
text(h100.XData(end),h100.YData(end),'100% increase') % 1% increase
text(hp1.XData(end),hp1.YData(end),'0.1% increase') % 1% increase
text(hnp1.XData(end),hnp1.YData(end),'0.1% decrease') % 1% increase
text(hn100.XData(end),hn100.YData(end),'100% decrease') % 1% increase

xlabel('$\dot\varepsilon$ $(a^{-1})$','interpreter','latex')
ylabel('$\Delta \dot\varepsilon$ $(a^{-1})$','interpreter','latex')
ytick=get(gca,'YTick');
yticklab=get(gca,'YTickLabel');
set(gca,'YTick',ytick(1:2:end),'YTickLabel',{yticklab{1:2:end}})
% }}}


%% 
figure(6);clf;
subplot(1,2,1);hold on;
scatter(le0(d<0),ld(d<0))
xlim([min(le0) max(le0)])
ylim([min(ld) max(ld)])
subplot(1,2,2);hold on;
scatter(le0(d>0),ld(d>0))
xlim([min(le0) max(le0)])
ylim([min(ld) max(ld)])

%plot(le0,ld,'.')
%plot(le0,abs(adresults.gradient),'.')
plot(le0,log10(abs(adresults.gradient)),'.')
plot(ld,log10(abs(adresults.gradient)),'.')
plot3(le0,ld,loggradn,'.')
figure(6);clf;hold on;
loggradn=log10(abs(gradn));
loggradn(gradn==0)=eps(max(gradn));
caxmin = (nanmean(loggradn) - 2*nanstd(loggradn));
caxmax = (nanmean(loggradn) + 3*nanstd(loggradn));
plotmodel(md,'data',loggradn,'caxis',[caxmin caxmax],'figure',6);

%
%
%% get eff stress
%ind=1; % from the first time step
%vx(vx==0)=NaN;
%vy(vy==0)=NaN;
%
%t0=md.results.TransientSolution(j).time;
%title(['vel change ' num2str(t0-5) '--' num2str(t0)])
%ctick=[-1000 -100 -10 -1 0 1 10 100 1000];
%logctick=sign(ctick).*log(1+abs(ctick)/10^C);
%h=get(gca,'colorbar');
%set(h,'Ticks',logctick,'TickLabels',ctick)
%h.Label.String='\Delta Vel (m/yr)';
%plotmodel(md,'data',effstrainrate,'caxis',[0 0.7E-8],'figure',2,'colormap','parula')
%
%n=md.materials.rheology_n; % Glen flow law exponent (non. dim.)
%B=md.materials.rheology_B; % rate factor (Pa s^(1/n) )
%effstress=B.*effstrainrate.^(1./n); % effective stress (Pa)
%
%% interpolate onto grid
%index=md.mesh.elements; % vertex indices of the mesh elements
%xlist=md.mesh.x(index); % x of vertecies of mesh elements
%ylist=md.mesh.y(index); % y of vertecies of mesh elements
%xc=mean(xlist,2); % element centroid in x
%yc=mean(ylist,2); % element centroid in y
%
%nx=1000; % number of points in x dim
%aspectratio=range(yc)/range(xc); % ny/nx
%xq=linspace(min(xc),max(xc),nx);
%yq=linspace(min(yc),max(yc),round(nx*aspectratio));
%[Xq,Yq]=meshgrid(xq,yq);
%
%F = scatteredInterpolant(xc,yc,effstress);
%effstressq=F(Xq,Yq);
%% }}}
%
%% mask ocean and outside of domain
%effstressq(effstressq<=0)=nan; % mask ocean
%dom=expread('Exp/ThwaitesPIGDomain2.exp'); % ISSM domain
%in=inpolygon(Xq,Yq,dom.x,dom.y);
%effstressq(~in)=nan;
%
%
%%plotmodel(md,'data',(effstress*1E-3)) % plot eff stress (kPa)
%thold=prctile(effstressq(:),98);
%thold=200E3;
%E=effstressq>thold;
%boundaries=bwboundaries(E,'TraceStyle','pixeledge');
%
%figure(1); clf;
%imagesc('XData',xq,'YData',yq,'CData',effstressq*1E-3,'alphadata',~isnan(effstressq)) % plot eff stress (kPa)
%%imshow(E) % plot eff stress (kPa)
%cb=colorbar;
%cb.Label.String='Effective Stress (kPa)';
%% threshold boundaries
%hold on;
%dx=mean(diff(xq));
%dy=mean(diff(yq));
%
%%for k=1:length(boundaries)
%%   b = boundaries{k};
%%	bx=(b(:,2)-1)/(length(xq)-1)*range(xq)+xq(1);
%%	by=(b(:,1)-1)/(length(yq)-1)*range(yq)+yq(1);
%%	%bx=xq(b(:,2));
%%	%by=yq(b(:,1));
%%   plot(bx,by,'r','LineWidth',0.5);
%%end
%% axis setup
%clim([0 350])
%set(gca,'YDir','normal')
%axis equal tight off
%
%xl=xlim;
%plot([xl(1),xl(1)+100E3],[-8 -8]*1E5,'k-','LineWidth',4)
%text(xl(1),-7.85E5,'100 km')
%
%
%%figure(2); clf; hold on;
%%plot(sort(effstressq(:)));
%%yline(thold);

