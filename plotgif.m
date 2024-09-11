for i=1:2
	if i==1
		md=loadmodel('Models/Amundsen_n3_csiro85_TransientRun.mat');
		fname='velocitychange3.gif';
	else
		md=loadmodel('Models/Amundsen_n4_csiro85_TransientRun.mat');
		fname='velocitychange4.gif';
	end

	j0=7
	for j=(j0+1):length(md.results.TransientSolution)
		dvel=md.results.TransientSolution(j).Vel-md.results.TransientSolution(j-1).Vel;
		C=-1;
		logdvel=sign(dvel).*log(1+abs(dvel)/10^C);
		plotmodel(md,'data',logdvel);
		caxis([-log(1+abs(1000)/10^C) log(1+abs(1000)/10^C)]);
		colormap('bluewhitered');
		t0=md.results.TransientSolution(j).time;
		title(['vel change ' num2str(t0-1) '--' num2str(t0)])
		ctick=[-1000 -100 -10 -1 0 1 10 100 1000];
		logctick=sign(ctick).*log(1+abs(ctick)/10^C);
		h=get(gca,'colorbar');
		set(h,'Ticks',logctick,'TickLabels',ctick)
		h.Label.String='\Delta Vel (m/yr)';
		axis off
		exportgraphics(gcf,fname,'Append',true);
	end
end
