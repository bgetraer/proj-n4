% Map of grounding line evolution 
%  subplot 1: average grounding line position of the n=3 models every 50 years
%  subplot 2: average grounding line position of the n=4 models every 50 years

readfrommodel=0;
readfromfile=0;
% read in data {{{
datafname='./oceanlevelsetresults.mat';
if readfrommodel
	t={};oceanlevelsets={}; % initialize blank cell arrays
	ismip_models={'miroc','noresm','hadgem','csiro','ipsl','ccsm'};
	n3cases=append('n3_',ismip_models,'85');
	n4cases=append('n4_',ismip_models,'85');
	cases={n3cases{:} n4cases{:}};
   fname_prefix='./Models/Amundsen_';
   fname_trans_suffix='_TransientRun';
   fname_pickup_suffix='_PickupTransient';
   for i=1:length(cases)
      % check if pickup transient exists
      if exist([fname_prefix cases{i} fname_pickup_suffix '.mat'])
         fname=[fname_prefix cases{i} fname_pickup_suffix '.mat'];
      else
          fname=[fname_prefix cases{i} fname_trans_suffix '.mat'];
      end

      disp(['   Loading transient solutions from ' fname]);
      md=loadmodel(fname); % Transient run for n=i
      t{i}=cell2mat({md.results.TransientSolution.time})';
      oceanlevelsets{i}=cell2mat({md.results.TransientSolution.MaskOceanLevelset})';
   end
	x = md.mesh.x;
	y = md.mesh.y;
	save(datafname,'t','oceanlevelsets','x','y','cases');
elseif readfromfile
	disp(['   Loading data from ' datafname]);
   load(datafname);
end% }}}

% interpolation grid in space
xq = linspace(min(x),max(x),500);
yq = linspace(min(y),max(y),500);
[Xq,Yq]= meshgrid(xq,yq);

% get bedmachine bed topography
%bed = interpBedmachineAntarctica(Xq,Yq,'bed');

% interpolate in time
tq = flip([2300:-50:2013]);
Vq = zeros(length(tq),length(x),length(oceanlevelsets));
for i = 1:length(oceanlevelsets)
	Vq(:,:,i) = interp1(t{i},oceanlevelsets{i},tq); % these are the level sets at our desired time steps	
end

% average level set at timestep
meanVq3 = mean(Vq(:,:,1:6),3); % mean for n=3 models
meanVq4 = mean(Vq(:,:,7:end),3); % mean for n=4 models

% for each timestep, plot the zero levelset contour (grounding line)
cmap = flip(brewermap(length(tq),'rdpu'));
clear hc
legstr = cell(size(tq));
figure(1);clf;
set(gcf,'position',[0 0 900 700]);
s1=subplot(1,2,1);hold on;
imagesc(bed,'Xdata',xq,'Ydata',yq);
s2=subplot(1,2,2);hold on;
imagesc(xq,yq,bed);

for i = 1:length(tq)
	% counting backwards to plot first contour on top
	j = length(tq)-i+1; 
	legstr{j} = num2str(tq(j)); 

	subplot(1,2,1);hold on;
	% interpolate levelset spatially onto regular grid
	W = meanVq3(j,:)'; % this scattered levelset
	Wq = griddata(x,y,W,Xq,Yq); % this gridded levelset
	levels = [0,0]; %zero levelset contour
	[~,hc(j)] = contour(xq,yq,Wq,levels,'color',cmap(j,:),'linewidth',1.5); % plot the contour
	subplot(1,2,2);hold on;
	% interpolate levelset spatially onto regular grid
   W = meanVq4(j,:)'; % this scattered levelset
   Wq = griddata(x,y,W,Xq,Yq); % this gridded levelset
   levels = [0,0]; %zero levelset contour
   contour(xq,yq,Wq,levels,'color',cmap(j,:),'linewidth',1.5); % plot the contour
end
set([s1,s2],'fontsize',12)
set(s1,'yticklabel',get(s1,'ytick')*1E-3,'xticklabel',get(s1,'xtick')*1E-3)
set(s2,'yticklabel',{},'xticklabel',get(s2,'xtick')*1E-3)
xlabel([s1,s2],'x (km)');
ylabel(s1,'y (km)');
set([s1,s2],'layer','top');
subplot(1,2,1);axis equal;
subplot(1,2,2);axis equal;
title(s1,'average of n=3 models');
title(s2,'average of n=4 models');
lg=legend(s1,hc,legstr,'Location','n','NumColumns',3);
title(lg,'grounding line');
cb=colorbar(s2,'location','north','color','k','AxisLocation','in');
title(cb,'bed elevation (m)')
cb.Position=[0.5819    0.74    0.3119    0.0305];

disp('printing figure');
exportgraphics(gcf,'figures/avg_gl.eps','BackgroundColor','none','Resolution',600,'ContentType','vector')
