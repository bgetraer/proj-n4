% Pot models of Amundsen Sea Embayment with n=3 or n=4 for comparison
% Benjamin Getraer
% Last Edited: 1/10/2024

%constants
rhoGt = 917*1e-12; % density of ice in Gt/m^3
gt2mmslr = 1/361.8;	% 361.8 Gt of ice will raise global sea levels by ~1 mm

ismip_models={'miroc','noresm','hadgem','csiro','ipsl','ccsm'};
n3cases=append('n3_',ismip_models,'85');
n4cases=append('n4_',ismip_models,'85');
cases={n3cases{:} n4cases{:}};

str = {'MIROC-ESM-CHEM','NorESM1-M','HadGEM2-ES','CSIRO-Mk3-6-0','IPSL-CM5A-MR','CCSM4'};
%s3=append(ismip_models,' (n=3)');
%s4=append(ismip_models,' (n=4)');
%s={s3{:} s4{:}};

readfrommodel=0;
% read in data {{{
datafname='./nThwaites.mat';
t={};v={};gt={}; % initialize blank cell arrays
if readfrommodel==1
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
		t{i}=cell2mat({md.results.TransientSolution(:).time})';
%		if any(i==[1:4])
%			t{i}=t{i}+2013;
%		end
		v{i}=cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation})';
		gt{i}=(v{i}-v{i}(1))*rhoGt;
%		A=[t{i},t{i}.^2,t{i}.^3];
%		m{i}=(A'*A)\A'*gt{i};
	end
	save(datafname,'t','v','gt','s');
else
	load(datafname);
end% }}}

fignum=5;
nmodels=length(ismip_models); % number of ISMIP models
fontname='Noto Mono';
switch fignum
	case 1 
		% FIG 1: N3 and N4 through 2100 {{{
		% color scheme ------------------------ {{{
		colors=[... 
			0.4000    0.7608    0.6471; ... % blue green 
			0.9882    0.5529    0.3843; ... % coral
			0.2196    0.4235    0.6902; ... % dk blue 
			1.0000    0.8510    0.1843; ... % yellow
			0.9412    0.0078    0.4980; ... % hot pink
			0.7451    0.6824    0.8314; ... % lavender
			]; % }}}
		% initalize figure and axes ---------------- {{{
		figure(1);clf;
		ax=axes;hold on;
		yyaxis left
		% }}}
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:20:2100]; % user defined grid X [start:spaces:end]
		gy_major=[-18:2:0]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+10; % user defined grid X [start:spaces:end]
		gy_minor=gy_major+1E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{
		%[~, ind]=sort(abs(sum(cell2mat(gt(1:nmodels))))); % sort models
		ind=1:6;
		for i=1:nmodels
			j=ind(i);
			h{i}=plot(t{j}(t{j}<=2100),gt{j}(t{j}<=2100),'-','color',colors(i,:),'linewidth',2);
			str{i}=s_leg{j};
		end
		for i=1:nmodels
			j=ind(i)+6;
			h{i+6}=plot(t{j}(t{j}<=2100),gt{j}(t{j}<=2100),'--','color',colors(i,:),'linewidth',2);
		end
		% }}}
		% axis options --------------------------- {{{

		xlim([2013 2100]);
		ylimleft=[-17E3 0];
		ylim(ylimleft);
		xlabel('time (yr)');
		ylabel('Change in mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0],'fontsize',12);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):5:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		ylabel('Sea level rise equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}
		% legend --------------------------------- {{{
		dh3=plot([1 2],[1,2],'k-','linewidth',2); % dummy handle for n=3 lines
		dh4=plot([1 2],[1,2],'k--','linewidth',2); % dummy handle for n=4 lines
		leg_ds3=['n=3{          }']; % dummy legend text for n=3 lines
		leg_ds4=['n=4{          }']; % dummy legend text for n=4 lines
		% first legend
		for i=1:length(str)
			spc_buff{i}=repmat(' ',1,7-length(str{i}));
		end
		leg_s=strcat(str,spc_buff,'rcp8.5');
		fontname='Noto Mono';
		leg1=legend([h{1:nmodels}],{leg_s{1:nmodels}},'location','sw','fontname',fontname);
		% axes for the second legend
		axd = axes('position',get(gca,'position'),'visible','off','fontsize',12);
		leg2=legend(axd,[dh3,dh4],{leg_ds3,leg_ds4},'location','sw','fontname',fontname);
		leg2.Position=leg2.Position+[0,leg1.Position(4)+0.01,leg1.Position(3)-leg2.Position(3),0];
		%leg2.Interpreter='latex';
		% }}}
		set(gcf,'color',[1 1 1]);
		%export_fig('N3N42100.png')
		clear h
		% }}}
	case 2 
		% FIG 2: Spread of models {{{
		% interpolate values onto single time vector
		t0=2013;
		tend=2300;
		t_q=t0:tend; % query time
		gt_q=zeros(length(nmodels*2),length(t_q));
		for i=1:nmodels*2
			gt_q(i,:)=interp1(t{i},gt{i},t_q);
			gt_q(i,1)=0;
		end
		gt_q3=gt_q(1:6,:);
		gt_q4=gt_q(7:12,:);

		diffn4=abs(gt_q4-gt_q3);

		figure(2); clf;hold on;
		set(gcf,'position',[0 0 900 600],'color',[1 1 1]);
		% inset highlight --------------------------- {{{
		bgcolor=[1 1 0.85]; % yellow background for inset
		bgcolor=[237 237 238]/250; % gray background for inset
		hr=rectangle('position',[0,0,2100,10E3],'facecolor',bgcolor,'edgecolor','none');
		% }}}
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:50:2300]; % user defined grid X [start:spaces:end]
		gy_major=[0:10:100]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+25; % user defined grid X [start:spaces:end]
		%gy_minor=gy_major+1E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			%	plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{
		%for i=1:nmodels
		%   j=ind(i);
		%   h{i}=plot(t_q,diffn4(j,:),'-','color',colors(i,:),'linewidth',2);
		%   str{i}=ismip_models{j};
		%end
		hrng=plot(t_q,range(gt_q3),'r-','linewidth',2);
		hstd=plot(t_q,std(gt_q3),'r--','linewidth',2);

		hmax=plot(t_q,max(diffn4),'k-','linewidth',2);
		hmean=plot(t_q,mean(diffn4),'k--','linewidth',2);

		% other measures of dispersion
		%plot(t_q,meanabsdiff(gt_q3),'b:','linewidth',2); % mean absolute difference
		%plot(t_q,mad(gt_q3),'b-','linewidth',2); % mean absolute deviation
		%plot(t_q,mad(gt_q3,1),'b--','linewidth',2); % median absolute deviation
		%leg_s_add={'MD(ISMIP_{n3})','mad(ISMIP_{n3})','mad1(ISMIP_{n3})'};
		% }}}
		% axis options --------------------------- {{{
		set(gca,'fontsize',12);
		xlim([2013 2300]);
		ylimleft=[0 100E3];
		ylim(ylimleft);
		xlabel('time (yr)');
		ylabel('Mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0],'fontsize',12);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):20:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		ylabel('Sea level equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}
		% legend --------------------------------- {{{
		legend([hrng,hstd,hmax,hmean],{'range of n=3 models','std of n=3 models',...
			'max abs diff between n=3 and n=4','mean abs diff between n=3 and n=4'},...
			'location','nw','fontname',fontname);
		% }}}
		% subaxes {{{
		ax2= axes('position',[0.175 0.31 0.35 0.4],'color',bgcolor);
		set(ax2,'position',[0.175 0.33 0.35 0.4]);
		hold on;
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:20:2100]; % user defined grid X [start:spaces:end]
		gy_major=[0:2:10]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+10; % user defined grid X [start:spaces:end]
		%gy_minor=gy_major+1E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			%	plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{
		hrng=plot(t_q(t_q<=2100),range(gt_q3(:,t_q<=2100)),'r-','linewidth',2);
		hstd=plot(t_q(t_q<=2100),std(gt_q3(:,t_q<=2100)),'r--','linewidth',2);

		hmax=plot(t_q(t_q<=2100),max(diffn4(:,t_q<=2100)),'k-','linewidth',2);
		hmean=plot(t_q(t_q<=2100),mean(diffn4(:,t_q<=2100)),'k--','linewidth',2);
		% }}}
		% axis options --------------------------- {{{
		set(gca,'fontsize',11);
		xlim([2013 2100]);
		ylimleft=[0 10E3];
		ylim(ylimleft);
		%xlabel('time (yr)');
		%ylabel('Spread in \Delta mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0]);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):5:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		%ylabel('Sea level rise equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}

		% }}}
		%export_fig('spread_n3n4.eps')
		% }}}
	case 3 
		% FIG 3: N3 and N4 through 2300 {{{
		% color scheme ------------------------ {{{
		colors=[... 
			0.4000    0.7608    0.6471; ... % blue green 
			0.9882    0.5529    0.3843; ... % coral
			0.2196    0.4235    0.6902; ... % dk blue 
			1.0000    0.8510    0.1843; ... % yellow
			0.9412    0.0078    0.4980; ... % hot pink
			0.7451    0.6824    0.8314; ... % lavender
			]; % }}}
		% initalize figure and axes ---------------- {{{
		figure(3);clf;
		set(gcf,'position',[0 0 900 600],'color',[1 1 1]);
		ax=axes;hold on;
		yyaxis left
		% }}}
		% inset highlight --------------------------- {{{
		bgcolor=[1 1 0.85]; % yellow background for inset
		bgcolor=[237 237 238]/250; % gray background for inset
		hr=rectangle('position',[0,-17.5E3,2100,17.5E3],'facecolor',bgcolor,'edgecolor','none');
		% }}}
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:50:2300]; % user defined grid X [start:spaces:end]
		gy_major=[-300:50:0]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+25; % user defined grid X [start:spaces:end]
		gy_minor=gy_major+25E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{
		for i=1:nmodels
			h{i}=plot(t{i}(t{i}<=2300),gt{i}(t{i}<=2300),'-','color',colors(i,:),'linewidth',2);
		end
		for i=1:nmodels
			h{i+6}=plot(t{i+6}(t{i+6}<=2300),gt{i+6}(t{i+6}<=2300),'--','color',colors(i,:),'linewidth',2);
		end
		% }}}
		% axis options --------------------------- {{{

		xlim([2013 2300]);
		ylimleft=[-275E3 0];
		ylim(ylimleft);
		xlabel('time (yr)');
		ylabel('\Delta mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0],'fontsize',12);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):100:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		ylabel('Sea level equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}
		% legend --------------------------------- {{{
		dh3=plot([1 2],[1,2],'k-','linewidth',2); % dummy handle for n=3 lines
		dh4=plot([1 2],[1,2],'k--','linewidth',2); % dummy handle for n=4 lines
		leg_ds3=['n=3']; % dummy legend text for n=3 lines
		leg_ds4=['n=4']; % dummy legend text for n=4 lines
		% first legend
		for i=1:length(str)
			spc_buff{i}=repmat(' ',1,7-length(str{i}));
		end
		leg_s=strcat(str,spc_buff);
		leg1=legend([h{1:nmodels}],{leg_s{1:nmodels}},'location','sw','fontname',fontname,'numcolumns',2);
		% axes for the second legend
		axd = axes('position',get(gca,'position'),'visible','off','fontsize',12);
		leg2=legend(axd,[dh3,dh4],{leg_ds3,leg_ds4},'location','sw','fontname',fontname);
		leg2.Position=leg2.Position+[leg1.Position(3)+0.01,0,0,leg1.Position(4)-leg2.Position(4)];
		% }}}
		% subaxes {{{
		ax2= axes('position',[0.175 0.31 0.36 0.42],'color',bgcolor);
		hold on;
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:20:2100]; % user defined grid X [start:spaces:end]
		gy_major=[-20:5:0]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+10; % user defined grid X [start:spaces:end]
		gy_minor=gy_major+2.5E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{
		%[~, ind]=sort(abs(sum(cell2mat(gt(1:nmodels))))); % sort models
		ind=1:6;
		for i=1:nmodels
			j=ind(i);
			h{i}=plot(t{j}(t{j}<=2100),gt{j}(t{j}<=2100),'-','color',colors(i,:),'linewidth',2);
			str{i}=ismip_models{j};
		end
		for i=1:nmodels
			j=ind(i)+6;
			h{i+6}=plot(t{j}(t{j}<=2100),gt{j}(t{j}<=2100),'--','color',colors(i,:),'linewidth',2);
			%str{i+6}=ismip_models{j};
		end
		% }}}
		% axis options --------------------------- {{{
		set(gca,'fontsize',11);
		xlim([2013 2100]);
		ylimleft=[-17.5E3 0];
		ylim(ylimleft);
		%xlabel('time (yr)');
		%ylabel('Spread in \Delta mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0]);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):5:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		%ylabel('Sea level rise equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}

		% }}}
		set(gcf,'color',[1 1 1]);
		%export_fig('N3N42300.eps')
		clear h

		% }}}
	case 4 
		% FIG 4: AD of n wrt mass above flotation {{{
		figure(4); clf;
		subplot(1,3,1); hold on;
		subplot(1,3,2); hold on;
		subplot(1,3,3); hold on;
		% }}}
	case 5 
		% FIG 5: mean N3 and N4 through 2300 {{{

		% interpolate values onto single time vector
		t0=2013;
		tend=2300;
		t_q=t0:tend; % query time
		gt_q=zeros(length(nmodels*2),length(t_q));
		for i=1:nmodels*2
			gt_q(i,:)=interp1(t{i},gt{i},t_q);
			gt_q(i,1)=0;
		end
		gt_q3=gt_q(1:6,:);
		gt_q4=gt_q(7:12,:);
		% initalize figure and axes ---------------- {{{
		figure(3);clf;
		set(gcf,'position',[0 0 900 600],'color',[1 1 1]);
		ax=axes;hold on;
		yyaxis left
		% }}}
		% inset highlight --------------------------- {{{
		bgcolor=[1 1 0.85]; % yellow background for inset
		bgcolor=[237 237 238]/250; % gray background for inset
		hr=rectangle('position',[0,-17.5E3,2100,17.5E3],'facecolor',bgcolor,'edgecolor','none');
		% }}}
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:50:2300]; % user defined grid X [start:spaces:end]
		gy_major=[-300:50:0]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+25; % user defined grid X [start:spaces:end]
		gy_minor=gy_major+25E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{		
			% shade min/max
			%h_range3 = patch([t_q,flip(t_q)],[min(gt_q3),flip(max(gt_q3))],'r','facealpha',0.25,'edgecolor','none');
			%h_range4 = patch([t_q,flip(t_q)],[min(gt_q4),flip(max(gt_q4))],'k','facealpha',0.25,'edgecolor','none');

			% plot mean lines
			h_mean3 = plot(t_q,mean(gt_q3),'-','color','r','linewidth',2);
			h_mean4 = plot(t_q,mean(gt_q4),'-','color','k','linewidth',2);

			% plot min/max
			h_range3=plot(t_q,min(gt_q3),'--','color','r','linewidth',2);
			plot(t_q,max(gt_q3),'--','color','r','linewidth',2);
			h_range4=plot(t_q,min(gt_q4),'--','color','k','linewidth',2);
			plot(t_q,max(gt_q4),'--','color','k','linewidth',2);

			% plot std
%			h{i}=plot(t_q,mean(gt_q3)-std(gt_q3),'--','color','r','linewidth',2);
%         h{i}=plot(t_q,mean(gt_q3)+std(gt_q3),'--','color','r','linewidth',2);
%         h{i}=plot(t_q,mean(gt_q4)-std(gt_q4),'--','color','k','linewidth',2);
%         h{i}=plot(t_q,mean(gt_q4)+std(gt_q4),'--','color','k','linewidth',2);
		% }}}
		% axis options --------------------------- {{{

		xlim([2013 2300]);
		ylimleft=[-275E3 0];
		ylim(ylimleft);
		xlabel('time (yr)');
		ylabel('\Delta mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0],'fontsize',12);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):100:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		ylabel('Sea level equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}
		% legend --------------------------------- {{{
		% first legend
		leg1=legend([h_mean3,h_mean4,h_range3,h_range4],{'mean and','mean and','range of n=3 models','range of n=4 models'},'location','sw','fontname',fontname,'numcolumns',2);
		% }}}
		% subaxes {{{
		ax2= axes('position',[0.175 0.31 0.36 0.42],'color',bgcolor);
		hold on;
		% gridlines --------------------------- {{{
		hold on
		gx_major=[2000:20:2100]; % user defined grid X [start:spaces:end]
		gy_major=[-20:5:0]*1E3; % user defined grid Y [start:spaces:end]
		gx_minor=gx_major+10; % user defined grid X [start:spaces:end]
		gy_minor=gy_major+2.5E3; % user defined grid Y [start:spaces:end]
		for i=1:length(gx_major)
			plot([gx_major(i) gx_major(i)],[gy_major(1) gy_major(end)],'-','color',[1,1,1]*0.8) %y grid lines
			plot([gx_minor(i) gx_minor(i)],[gy_major(1) gy_major(end)],':','color',[1,1,1]*0.8) %y grid lines
		end
		for i=1:length(gy_major)
			plot([gx_major(1) gx_major(end)],[gy_major(i) gy_major(i)],'-','color',[1,1,1]*0.8) %x grid lines
			plot([gx_major(1) gx_major(end)],[gy_minor(i) gy_minor(i)],':','color',[1,1,1]*0.8) %x grid lines
		end
		% }}}
		% plot data ------------------------------ {{{
			% shade min/max
			%patch([t_q,flip(t_q)],[min(gt_q3),flip(max(gt_q3))],'r','facealpha',0.25,'edgecolor','none');
			%patch([t_q,flip(t_q)],[min(gt_q4),flip(max(gt_q4))],'k','facealpha',0.25,'edgecolor','none');

			% plot mean lines
			plot(t_q,mean(gt_q3),'-','color','r','linewidth',2);
			plot(t_q,mean(gt_q4),'-','color','k','linewidth',2);
			% plot min/max
			plot(t_q,min(gt_q3),'--','color','r','linewidth',2);
			plot(t_q,max(gt_q3),'--','color','r','linewidth',2);
			plot(t_q,min(gt_q4),'--','color','k','linewidth',2);
			plot(t_q,max(gt_q4),'--','color','k','linewidth',2);
		% }}}
		% axis options --------------------------- {{{
		set(gca,'fontsize',11);
		xlim([2013 2100]);
		ylimleft=[-17.5E3 0];
		ylim(ylimleft);
		%xlabel('time (yr)');
		%ylabel('Spread in \Delta mass above flotation (10^3 Gt)');
		set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0]);
		yyaxis right
		ylimright=ylimleft.*gt2mmslr;
		ytickright=ylimleft(1):5:ylimleft(2);
		set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
		set(gca,'YTickLabel',abs(get(gca,'YTick')));
		%ylabel('Sea level rise equivalent (mm)');
		set(gca, 'Layer', 'top');
		% }}}

		% }}}
		set(gcf,'color',[1 1 1]);
		%export_fig('MeanN3N4.eps')
		clear h

		% }}}
end 

% How much more ice if n=4? {{{
for i=1:6
	n3_2100(i)=interp1(t{i},gt{i},2100);
	n4_2100(i)=interp1(t{i+6},gt{i+6},2100);

	n3_2300(i)=interp1(t{i},gt{i},2300);
   n4_2300(i)=interp1(t{i+6},gt{i+6},2300);
end
m2100 = mean((n4_2100-n3_2100)./n3_2100);
s2100 = std((n4_2100-n3_2100)./n3_2100);
m2300 = mean((n4_2300-n3_2300)./n3_2300);
s2300 = std((n4_2300-n3_2300)./n3_2300);
disp('Mass loss by 2100:');
fprintf('   %0.2f (%0.2f) times more ice if n=4 \n',m2100,s2100);
disp('Mass loss by 2300:');
fprintf('   %0.2f (%0.2f) times more ice if n=4 \n',m2300,s2300);
% }}}
%% Figure 2 {{{
%figure(2);clf;hold on;
%for i=n
%	gt0=interp1(t{3},gt{3},t{i});
%   h{i}=plot(t{i},(gt{i}-gt0),'color',c(i,:),'linewidth',2);
%end
%%h{3}=plot(xlim,[0, 0],'color',c{3},'linewidth',2);
%legend([h{:}],sl,'location','sw');
%xlim([2013 2100]);
%xlabel('time (yr)');
%ylabel('Mass anomaly (10^3 Gt)');
%set(gca,'YTickLabel',get(gca,'YTick')*1E-3,'YColor',[0,0,0]);
%ylimleft=get(gca,'YLim');
%yyaxis right
%ylimright=ylimleft.*gt2mmslr;
%ytickright=ylimleft(1):50:ylimleft(2);
%set(gca,'YTick',ytickright,'Ylim',ylimright,'YColor',[0,0,0]);
%ylabel('SLR anomaly (mm)');
%title('Thwaites Glacier mass projection anomaly for different n')
%% }}}



%%read in ismip data {{{
%readfrommodel=0;
%datafname='./ismipThwaites.mat';
%if readfrommodel==1
%	transfile='./Models/Amundsen_n3_%s_TransientRun';
%	fname={sprintf(transfile,ismip{1}),sprintf(transfile,ismip{2})};
%	for i=1:length(fname)
%		disp(['   Loading transient solutions from ' fname{i}]);
%		md=loadmodel(fname{i}); % Transient run for n=i
%		t{i}=cell2mat({md.results.TransientSolution(:).time})';
%		v{i}=cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation})';
%		gt{i}=(v{i}-v{i}(1))*rhoGt;
%		A=[t{i},t{i}.^2,t{i}.^3];
%		m{i}=(A'*A)\A'*gt{i};
%		s{i}=sprintf('n=%i',i);
%	end
%	save(datafname,'t','v','gt','m','s');
%else
%	%load(datafname);
%end% }}}
%%load transientsolutions {{{ 
%if ~exist('md3') | ~exist('md4')
%	disp('   Loading transient solutions');
%	filestr = './Models/Amundsen_n%0.0d_TransientRun';
%	
%	% n=3 model 
%	md3 = loadmodel(sprintf(filestr,3)); % Transient run for n=3
%	t3=cell2mat({md3.results.TransientSolution(:).time})';
%	v3=cell2mat({md3.results.TransientSolution(:).IceVolumeAboveFloatation})';
%	% n=4 model
%	md4 = loadmodel(sprintf(filestr,4)); % Transient run for n=4
%	t4=cell2mat({md4.results.TransientSolution(:).time})';
%	v4=cell2mat({md4.results.TransientSolution(:).IceVolumeAboveFloatation})';
%end % }}}
%%plot initial velocities {{{
%%plotmodel(md3,'data',md3.initialization.vel-md4.initialization.vel,'figure',1);
%% }}}
%%plot change in volume {{{
%
%%constants
%mcube2gt = md3.materials.rho_ice*1e-12;
%gt2mmslr = 1/360;
%
%%volume to gigatons
%gt3=(v3-v3(1))*mcube2gt;
%gt4=(v4-v4(1))*mcube2gt;
%
%A3 = [ t3, t3.^2];
%m3 = (A3'*A3)\A3'*gt3;
%A4 = [ t4, t4.^2];
%m4 = (A4'*A4)\A4'*gt4;
%
%figure(2);clf;hold on;
%h3 = plot(t3,gt3,'.');
%h4 = plot(t4,gt4,'.');
%
%lg3 = sprintf('n=3: y\\approx%0.1ft^2+%3.0ft',m3([2,1]));
%lg4 = sprintf('n=4: y\\approx%0.1ft^2+%3.0ft',m4([2,1]));
%ylabel('mass (Gt)')
%xlabel('t (years)')
%
%plot(xlim,gt3(end)*[1,1],'--k')
%plot(xlim,gt4(end)*[1,1],'--k')
%dgt = gt4(end)-gt3(end);
%text(20,gt4(end) - 1/2*(dgt),...
%	sprintf('\\Delta_{100}=%0.0f Gt \\equiv %0.0fmm slr',[dgt,dgt*gt2mmslr]))
%legend([h3,h4],lg3,lg4);
%
%%figure(3);clf;hold on;
%%%hdiff = plot(t,gt4-gt3,'.');
%%hmdiff = plot(t,A4*m4-A3*m3));
%%
%%lgdiff = 'n=4 - n=3';
%%lgmdiff = sprintf('y\\approx%0.1ft^2+%3.0ft',m4([2,1])-m3([2,1]));
%%legend([hdiff,hmdiff],lgdiff,lgmdiff);
%%ylabel('mass (Gt)');
%%xlabel('t (years)');
%%title('Thwaites Glacier: difference between n=4 and n=3 models')
%%
%%dmap3 = md3.results.TransientSolution(end).Thickness - ...
%%md3.results.TransientSolution(1).Thickness;
%%dmap4 = md4.results.TransientSolution(end).Thickness - ...
%%md4.results.TransientSolution(1).Thickness;
%%plotmodel(md3,'data',dmap3,'figure',4)
%%plotmodel(md4,'data',dmap4,'figure',5)
%%plotmodel(md4,'data',dmap4-dmap3,'figure',6)
%
%% }}}

%for j=(1):length(md.results.TransientSolution)
%	dvel=md.results.TransientSolution(end).Vel-md.results.TransientSolution(j).Vel;
%
%	for j=(5+1):length(md.results.TransientSolution)
%		dvel=md.results.TransientSolution(j).Vel-md.results.TransientSolution(j-5).Vel;
%		C=-1;
%		logdvel=sign(dvel).*log(1+abs(dvel)/10^C);
%		plotmodel(md,'data',logdvel);
%		caxis([-log(1+abs(1000)/10^C) log(1+abs(1000)/10^C)]);
%		colormap('bluewhitered');
%		t0=md.results.TransientSolution(j).time;
%		title(['vel change ' num2str(t0-5) '--' num2str(t0)])
%		ctick=[-1000 -100 -10 -1 0 1 10 100 1000];
%		logctick=sign(ctick).*log(1+abs(ctick)/10^C);
%		h=get(gca,'colorbar');
%		set(h,'Ticks',logctick,'TickLabels',ctick)
%		h.Label.String='\Delta Vel (m/yr)';
%		axis off
%		%exportgraphics(gcf,'velocitychange.gif','Append',true);
%	end
%end
