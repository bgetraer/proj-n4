% Generate models of Amundsen Sea Embayment with variable n and different ISMIP6 forcing models for comparison
%    originating from a course project for Dartmouth EARS107 (Fall 2022).
%
%    A single model is initialized using the n=4 assumptions to invert for B element-wise. After inverting for B and C,
% the stress balance is solved to generate the initial model velocity field. We then diverge the models being tested, 
% converting the B values for n=4 such that the initial model velocity and viscocity fields are EXACTLY the same between 
% the different models. 
%    By default, atmospheric forcing is taken from RACMO and melting rates from Rignot et al. If ``useISMIP6'' is turned on,
% these forcings are replaced by interpolations of the chosen ISMIP6 model. 

% SEE ALSO interpISMIP6AntarcticaOcn,  interpISMIP6AntarcticaSMB
% Last Edited: 1/25/2024

steps = [3:5];

rheology_n=3;  % Glen's flow law exponent used for this experiment
useISMIP6=0;	% turn on ISMIP6 rcp 8.5 climate scenario? otherwise use RACMO
% define ISMIP model {{{
ISMIP6model.name='ipsl85';  % name of climate model
if useISMIP6
	switch ISMIP6model.name
		case 'miroc85'
			ISMIP6model.oceanmodelname='miroc-esm-chem_rcp8.5';
			ISMIP6model.atmosmodelname='miroc-esm-chem_rcp8.5';
		case 'noresm85'
			ISMIP6model.oceanmodelname='noresm1-m_rcp8.5';
			ISMIP6model.atmosmodelname='noresm1-m_rcp8.5';
		case 'ccsm85'
			ISMIP6model.oceanmodelname='ccsm4_rcp8.5';
			ISMIP6model.atmosmodelname='ccsm4_rcp8.5';
		case 'hadgem85'
			ISMIP6model.oceanmodelname='hadgem2-es_rcp8.5';
			ISMIP6model.atmosmodelname='HadGEM2-ES_rcp85';
		case 'csiro85'
			ISMIP6model.oceanmodelname='csiro-mk3-6-0_rcp8.5';
			ISMIP6model.atmosmodelname='CSIRO-Mk3-6-0_rcp85';
		case 'ipsl85'
			ISMIP6model.oceanmodelname='ipsl-cm5a-mr_rcp8.5';
			ISMIP6model.atmosmodelname='IPSL-CM5A-MR_rcp85';
		otherwise
			error('ISMIP6 model name not supported')
	end
	ISMIP6model.prefix=[ISMIP6model.name '_'];
else 
	ISMIP6model.prefix=''; % make blank for the file naming scheme
end
% }}}

% cluster parameters{{{
cluster=generic('name',oshostname(),'np',45); %for totten 45 ideal
%}}}
% model naming scheme{{{ 
prefix_initial = 'Amundsen_backwards_'; %ie 'Amundsen_Mesh.mat'
prefix_n = [prefix_initial 'n' num2str(rheology_n) '_']; %ie 'Amundsen_n4_TransientRun.mat'
prefix_trans = [prefix_n ISMIP6model.prefix];
%}}}

org=organizer('repository',['./Models'],'prefix',prefix_initial,'steps',steps); clear steps;
addpath([getenv('JPL_DIR') '/proj-morlighem/CODE/']);
addpath('./m/');

%Model initialization
if perform(org,'InversionB'),% {{{
   md=loadmodel('Models/Amundsen_Param.mat');

	% set n and A
	md.materials.rheology_n = 4.*ones(md.mesh.numberofelements,1);
	md.materials.rheology_B = 22631144 .* ones(md.mesh.numberofelements,1);

   % set M1QN3 package
   md.inversion=m1qn3inversion(md.inversion);

   % Set inversion data
   md.inversion.vx_obs=md.initialization.vx; % initialization was defined in last step (raw data, with NaN)
   md.inversion.vy_obs=md.initialization.vy; % initialization was defined in last step (raw data, with NaN)
   %Fill in blanks in velocity data
   pos=find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs));
   md.inversion.vx_obs(pos)=0;
   md.inversion.vy_obs(pos)=0;
   md.inversion.vel_obs=sqrt(md.inversion.vx_obs.^2+md.inversion.vy_obs.^2);
   md.initialization.vx(pos)=0;
   md.initialization.vy(pos)=0;
   md.initialization.vel(pos)=0;

   % Control general
   md.inversion.iscontrol=1;
   md.inversion.maxsteps=40;
   md.inversion.maxiter=40;
   md.inversion.dxmin=0.1;
   md.inversion.gttol=1.0e-6;
   md.inversion.incomplete_adjoint=0; % 0: non linear viscosity, 1: linear viscosity 04/29/2019 changed to non linear

   % Cost functions
   md.inversion.cost_functions=[101 103]; %removed 502 because we are working element-wise
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,length(md.inversion.cost_functions));
   md.inversion.cost_functions_coefficients(:,1)=73.0; %change this such that the cost functions 101 and 103 have the same contribution at the end of the inversion
   md.inversion.cost_functions_coefficients(:,2)=1; % always 1
   md.inversion.cost_functions_coefficients(pos,:)=0; % positions with NaN in the velocity data set

   % Force "mounds" as Dirichlet BC
   % This avoids problems during the stressbalance convergence solver
   incontour=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,expread('./Exp/Mound_region_1995_correction_rev00.exp'),'node',1);
   pos=find(incontour);
   md.stressbalance.spcvx(pos)=md.initialization.vx(pos);
   md.stressbalance.spcvy(pos)=md.initialization.vy(pos);
   md.stressbalance.spcvz(pos)=md.initialization.vz(pos);
	 % second region
   incontour=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,expread('./Exp/Mound_region2_1995_correction_rev00.exp'),'node',1);
   pos=find(incontour);
   md.stressbalance.spcvx(pos)=md.initialization.vx(pos);
   md.stressbalance.spcvy(pos)=md.initialization.vy(pos);
   md.stressbalance.spcvz(pos)=md.initialization.vz(pos);
   % third region (floating)
   incontour=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,expread('./Exp/Floating_region_1995_correction_rev00.exp'),'node',1);
   pos=find(incontour);
   md.mask.ice_levelset(pos)=1; % forcing no ice

   % Controls
   % setting initial guess for rheology B
   md.inversion.control_parameters={'MaterialsRheologyBbar'};
   % Benjy added 20% rescaling of min and max
   md.inversion.min_parameters=7958362.*ones(size(md.materials.rheology_B)); % from Seroussi et al, 2014
   md.inversion.max_parameters=141898341*ones(size(md.materials.rheology_B)); % from Seroussi et al, 2014

   % Additional parameters
   md.stressbalance.restol=0.0001; % 04/29/2019
   %md.stressbalance.reltol=0.01; % 04/29/2019
   %md.stressbalance.abstol=10; % 04/29/2019
   %md.stressbalance.maxiter=40; % 08/21/2019
   mds.stressbalance.maxiter=50; % 10/24/2019
   mds.stressbalance.reltol=NaN; % 11/05/2019
   mds.stressbalance.abstol=NaN; % 11/05/2019

   % Prepare to solve
   md.cluster=cluster;
   md.verbose=verbose('solution',false,'control',true);
   md.miscellaneous.name='inversion_B';
   mds=extract(md,md.mask.ocean_levelset<0);
   mds.friction.coefficient(:)=0; % make sure there is no basal friction

   % Solve
   mds.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
   mds.toolkits.DefaultAnalysis.ksp_max_it=500;
   mds.settings.solver_residue_threshold=1e-6;

   mds=solve(mds,'Stressbalance'); % only extracted model

   % Update model rheology_B accordingly
   md.materials.rheology_B(mds.mesh.extractedelements) = mds.results.StressbalanceSolution.MaterialsRheologyBbar;

   % PLOTTING
   figure(1);clf;hold on;
   pos = md.mask.ocean_levelset<0 & md.mask.ice_levelset<0;
   hist(md.materials.rheology_B(md.materials.rheology_B~=mode(md.materials.rheology_B)),100);
   ylim([0,3000]);
   ylabel('number of elements');
   xlabel('rheology_B');

   savemodel(org,md);
end
%}}}
if perform(org,'InversionC'),% {{{
	md=loadmodel(org,'InversionB'); % load the rescaled n=4 model

	% positions with NaN in the velocity data set get a weight of zero
	pos=find(isnan(md.inversion.vx_obs) | isnan(md.inversion.vy_obs));

	% Friction Law
	disp('   -- Setting up a Budd''s sliding law (with m=1, q=1)');
	md.friction.p=ones(md.mesh.numberofelements,1);
	md.friction.q=ones(md.mesh.numberofelements,1);
	md.inversion.cost_functions=[101 103 501];
	md.inversion.cost_functions_coefficients = ones(md.mesh.numberofvertices,3); % IMPORTANT: set as 0 again the vertices with no data (NaN)
	md.inversion.cost_functions_coefficients(:,1)=3.5e+2;% for budd change this such that the cost functions 101 and 103 have the same contribution
	md.inversion.cost_functions_coefficients(:,2)=1; % always 1
	md.inversion.cost_functions_coefficients(:,end)=1e-8;  % for budd
	md.inversion.cost_functions_coefficients(pos,:)=0; % positions with NaN in the velocity data set
	%Controls
	md.friction.coefficient=EstimateFric_Budd(md); % initial guess from Driving Stress (using r=1 and s=1)
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.min_parameters=0*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=3000*ones(md.mesh.numberofvertices,1);
	md.miscellaneous.name='inversion_drag_budd';

	% Control
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=80;
	md.inversion.maxiter=80;
	md.inversion.dxmin=0.01;
	md.inversion.gttol=1.0e-5;
	md.inversion.incomplete_adjoint=0; % 0: non linear viscosity, 1: linear viscosity 04/29/2019 changed to non linear

	% Keep friction 0 in the floating part
	pos=find(md.mask.ocean_levelset<0);
	md.friction.coefficient(pos)=0;

	% Additional parameters for stress balance
	md.stressbalance.restol=0.002; % 04/29/2019
	%md.stressbalance.reltol=0.01; % 04/29/2019
	%md.stressbalance.abstol=10; % 04/29/2019
	%md.stressbalance.maxiter=30; % 04/29/2019
	% updated?
	md.stressbalance.maxiter=50; % 10/24/2019
	md.stressbalance.reltol=NaN; % 11/05/2019
	md.stressbalance.abstol=NaN; % 11/05/2019

	% Set cluster
	md.cluster=cluster;
	md.verbose=verbose('solution',false,'control',true);

	% Define solver
	md.toolkits.DefaultAnalysis=bcgslbjacobioptions();% biconjugate gradient with block Jacobi preconditioner
	md.toolkits.DefaultAnalysis.ksp_max_it = 500;

	md=solve(md,'Stressbalance');

	% Update friction coefficient and velocity field
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
	md.initialization.vx= md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;
	md.initialization.vel=md.results.StressbalanceSolution.Vel;

	% Get the modeled velocity field that N3 and N4 models have in common	
	md.inversion.iscontrol=0;
	md.verbose.convergence = 1;
	md = solve(md,'sb');
	% Update velocity field
	md.initialization.vx= md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;
	md.initialization.vel=md.results.StressbalanceSolution.Vel;

	savemodel(org,md);
end%}}}
if perform(org,'ConvertN'), % {{{
	md=loadmodel(org,'InversionC');  

	% extract B from the existing model (calculated for n=4)
	B4 = md.materials.rheology_B;
	% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
	[strainrate] = strainrate_SSA(md,md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy,0);

	% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
	disp(['Converting B for n=' num2str(rheology_n)]);
	Bn = B4.*(strainrate.eff/md.constants.yts).^((rheology_n-1)/rheology_n - 3.0/4.0);

	% update with new parameters for n for all ice vertices
	ind=strainrate.eff>0; % index the ice, not the ocean, ignore zero strainrate
	md.materials.rheology_B(ind)=Bn(ind); % update B
	md.materials.rheology_n = rheology_n * ones(size(md.materials.rheology_n)); % set n

	% test if the  model immediately converges to the same modeled velocity solution
	testconvergence = 1;
	if testconvergence
		fprintf(['\n' ...
			'*********************************************************\n' ...		
			'TESTING CONVERGENCE: solution should converge in one step\n' ...
			'*********************************************************\n']);
		md.inversion.iscontrol=0;
		md.verbose.convergence = 1;
		mdtest = solve(md,'sb');
		clear mdtest;
	end	
	% Save the model
	org.prefix=prefix_n;
	savemodel(org,md);
end%}}}
if perform(org,'TransientPrep'),% {{{
	org.prefix=prefix_n;
	md=loadmodel(org,'ConvertN');

	% TRANSIENT PREP {{{
	%Set parameters
	md.inversion.iscontrol=0;
	md.settings.output_frequency = 20;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;
	md.timestepping.final_time=300;

	% We set the transient parameters
	md.transient.ismovingfront=0;
	md.transient.isthermal=0;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.isgroundingline=1;
	md.groundingline.migration = 'SubelementMigration';

	%Basal melt rate
	md.basalforcings=linearbasalforcings();
	md.basalforcings.deepwater_melting_rate=50.; % m/yr ice equivalent
	md.basalforcings.deepwater_elevation=-500;
	md.basalforcings.upperwater_melting_rate=0; % no melting for zb>=0
	md.basalforcings.upperwater_elevation=0; % sea level
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1); % no melting on grounded ice

	%Same as for the drag inversion
	md.stressbalance.restol=0.002;
	md.verbose.solution=1;
	md.verbose.convergence = 0;
	md.cluster = cluster;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};
	md.stressbalance.requested_outputs={'default'};
	% }}}
	% LOAD ISMIP6 MODEL IF REQUIRED {{{
	if useISMIP6
		% load forcings
		md.basalforcings=interpISMIP6AntarcticaOcn(md,ISMIP6model.oceanmodelname);
		md.smb=interpISMIP6AntarcticaSMB(md,ISMIP6model.atmosmodelname);
		% reload geothermal heat flux
		md.basalforcings.geothermalflux  = interpSeaRISE(md.mesh.x,md.mesh.y,'bheatflx_shapiro',-1);
		% set islocal to 0
		md.basalforcings.islocal=0;

		% adjust time correctly
		md.timestepping.start_time=2013; % nominal year for BedMachine
		md.timestepping.final_time=2100; % end year for ISMIP6 forcings 
	end
	% }}}

	org.prefix=prefix_trans;
	savemodel(org,md);
end%}}}
if perform(org,'TransientRun'),% {{{
	org.prefix=prefix_trans;
	md=loadmodel(org,'TransientPrep');

	% SOLVE
	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster = cluster;
	loadonly=0;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	savemodel(org,md);
end%}}}
if perform(org,'PickupTransient'), % {{{
	org.prefix=prefix_trans;
	md=loadmodel(org,'TransientRun');

	% retain previous transient solution
	lasttransientsolution=md.results.TransientSolution;

	% start time at last end
	start1=md.timestepping.start_time; % overall start
	end1=md.timestepping.final_time;   % first run end
	md.timestepping.start_time=end1;   % next run start
	md.timestepping.final_time=start1+300; % overall end
	disp(['restarting run at t=' num2str(md.timestepping.start_time)]);

	% reinitialize
	md.geometry.base = md.results.TransientSolution(end).Base;
	md.geometry.thickness = md.results.TransientSolution(end).Thickness;
	md.geometry.surface = md.results.TransientSolution(end).Surface;
	% reset the masks
	md.mask.ocean_levelset = md.results.TransientSolution(end).MaskOceanLevelset;
	% reset the initialization
	md.initialization.vx = md.results.TransientSolution(end).Vx;
	md.initialization.vy = md.results.TransientSolution(end).Vy;
	md.initialization.vel = md.results.TransientSolution(end).Vel;
	md.initialization.pressure = md.results.TransientSolution(end).Pressure;

	md.miscellaneous.name = [org.prefix org.steps(org.currentstep).string];
	md.cluster=cluster;
	md.cluster.interactive=0;
	md=solve(md,'tr','runtimename',0,'loadonly',1);

	% append new transient solution to the last one
	md.results.TransientSolution=[lasttransientsolution md.results.TransientSolution];
	% save model
	savemodel(org,md);
end % }}}

% OTHER TESTING 
if perform(org,'PerturbN_baseline'), % {{{
	md=loadmodel('Models/Amundsen_n3_ccsm85_TransientPrep');

	% run for 20 years
	md.settings.output_frequency=10;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;
	md.timestepping.final_time=15;

	% SOLVE
	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster = cluster;
	loadonly=0;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	savemodel(org,md);
end % }}}
if perform(org,'PerturbN_gradmin'), % {{{
	md=loadmodel('Models/Amundsen_n3_ccsm85_TransientPrep');
	gradN=getfield(load('Models/gradient.mat'),'gradient');

	% extract B from the existing model (calculated for n=3)
	B3 = md.materials.rheology_B;
	% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
	[strainrate] = strainrate_SSA(md,md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy,0);

	% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
	disp('Converting B for n=4');
	B4 = B3.*(strainrate.eff/md.constants.yts).^((4.0-1)/4.0 - 2.0/3.0);

	% update with new parameters for n for all ice vertices
	ind=gradN<prctile(gradN,1); % target the most negative values of gradient (lower than 1st percentile)
	md.materials.rheology_B(ind)=B4(ind); % update B only for the target elements
	md.materials.rheology_n(ind)=4; % set n=4 for the target elements

	% run for 20 years
	md.settings.output_frequency=10;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;
	md.timestepping.final_time=15;

	% SOLVE
	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster = cluster;
	loadonly=0;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	savemodel(org,md);

	t=cell2mat({md.results.TransientSolution.time});
	vaf=cell2mat({md.results.TransientSolution.IceVolumeAboveFloatation});
	plot(t,vaf);
end % }}}
if perform(org,'PerturbN_gradmax'), % {{{
	md=loadmodel('Models/Amundsen_n3_ccsm85_TransientPrep');
	gradN=getfield(load('Models/gradient.mat'),'gradient');

	% extract B from the existing model (calculated for n=3)
	B3 = md.materials.rheology_B;
	% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
	[strainrate] = strainrate_SSA(md,md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy,0);

	% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
	disp('Converting B for n=4');
	B4 = B3.*(strainrate.eff/md.constants.yts).^((4.0-1)/4.0 - 2.0/3.0);

	% update with new parameters for n for all ice vertices
	ind=gradN>prctile(gradN,99); % target the most positive values of gradient (higher than 99th percentile)
	md.materials.rheology_B(ind)=B4(ind); % update B only for the target elements
	md.materials.rheology_n(ind)=4; % set n=4 for the target elements

	% run for 20 years
	md.settings.output_frequency=10;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;
	md.timestepping.final_time=15;

	% SOLVE
	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster = cluster;
	loadonly=0;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	savemodel(org,md);

	t=cell2mat({md.results.TransientSolution.time});
	vaf=cell2mat({md.results.TransientSolution.IceVolumeAboveFloatation});
	plot(t,vaf);
end % }}}
if perform(org,'PerturbN_gradmin_1'), % {{{
	md=loadmodel('Models/Amundsen_n3_ccsm85_TransientPrep');
	gradN=getfield(load('Models/gradient.mat'),'gradient');

	% extract B from the existing model (calculated for n=3)
	B3 = md.materials.rheology_B;
	% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
	[strainrate] = strainrate_SSA(md,md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy,0);

	% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
	disp('Converting B for n=4');
	B4 = B3.*(strainrate.eff/md.constants.yts).^((4.0-1)/4.0 - 2.0/3.0);

	% update with new parameters for n for all ice vertices
	ind=find(gradN==min(gradN)); % target the element with the most negative gradient
	md.materials.rheology_B(ind)=B4(ind); % update B only for the target element
	md.materials.rheology_n(ind)=4; % set n=4 for the target element

	% run for 20 years
	md.settings.output_frequency=10;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;
	md.timestepping.final_time=15;

	% SOLVE
	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster = cluster;
	loadonly=0;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	savemodel(org,md);

	tmin1=cell2mat({md.results.TransientSolution.time});
	vafmin1=cell2mat({md.results.TransientSolution.IceVolumeAboveFloatation});
	plot(tmin1,vafmin1);
end % }}}
if perform(org,'PerturbN_gradmax_1'), % {{{
	md=loadmodel('Models/Amundsen_n3_ccsm85_TransientPrep');
	gradN=getfield(load('Models/gradient.mat'),'gradient');

	% extract B from the existing model (calculated for n=3)
	B3 = md.materials.rheology_B;
	% calculate effective strain rate from the MODELED velocities, do NOT average over the elements
	[strainrate] = strainrate_SSA(md,md.results.StressbalanceSolution.Vx, md.results.StressbalanceSolution.Vy,0);

	% calculate approximation for B, where mu(3) = mu(n), ie B3/eps_e^(2/3) = Bn/eps_e^((n-1)/n)
	disp('Converting B for n=4');
	B4 = B3.*(strainrate.eff/md.constants.yts).^((4.0-1)/4.0 - 2.0/3.0);

	% update with new parameters for n for all ice vertices
	ind=find(gradN==max(gradN)); % target the element with the most positive gradient
	md.materials.rheology_B(ind)=B4(ind); % update B only for the target element
	md.materials.rheology_n(ind)=4; % set n=4 for the target element

	% run for 20 years
	md.settings.output_frequency=10;
	md.timestepping=timesteppingadaptive();
	md.timestepping.time_step_max=0.05;
	md.timestepping.time_step_min=0.0005;
	md.timestepping.start_time=0;
	md.timestepping.final_time=15;

	% SOLVE
	md.miscellaneous.name = org.steps(org.currentstep).string;
	md.cluster = cluster;
	loadonly=0;
	md=solve(md,'tr','runtimename',0,'loadonly',loadonly);

	savemodel(org,md);

	tmax1=cell2mat({md.results.TransientSolution.time});
	vafmax1=cell2mat({md.results.TransientSolution.IceVolumeAboveFloatation});
	plot(tmax1,vafmax1);
end % }}}
if perform(org,'PlotPerturbN'), % {{{
	md=loadmodel('Models/Amundsen_PerturbN_baseline.mat');
	minmd=loadmodel('Models/Amundsen_PerturbN_gradmin.mat');
	maxmd=loadmodel('Models/Amundsen_PerturbN_gradmax.mat');
	min1md=loadmodel('Models/Amundsen_PerturbN_gradmin_1.mat');
	max1md=loadmodel('Models/Amundsen_PerturbN_gradmax_1.mat');

	t0=cell2mat({md.results.TransientSolution.time});
	vaf0=cell2mat({md.results.TransientSolution.IceVolumeAboveFloatation});
	tmin=cell2mat({minmd.results.TransientSolution.time});
	tmax=cell2mat({maxmd.results.TransientSolution.time});
	vafmin=cell2mat({minmd.results.TransientSolution.IceVolumeAboveFloatation});
	vafmax=cell2mat({maxmd.results.TransientSolution.IceVolumeAboveFloatation});
	tmin1=cell2mat({min1md.results.TransientSolution.time});
	vafmin1=cell2mat({min1md.results.TransientSolution.IceVolumeAboveFloatation});
	tmax1=cell2mat({max1md.results.TransientSolution.time});
	vafmax1=cell2mat({max1md.results.TransientSolution.IceVolumeAboveFloatation});

	figure(1); clf; hold on;

	%plot(t0,vaf0-vaf0(1),'linewidth',2,'k')
	vaf0=vaf0-vaf0(1);
	plot(tmin,vafmin-vafmin(1)-vaf0,'linewidth',2);
	plot(tmax,vafmax-vafmax(1)-vaf0,'linewidth',2);
	plot(tmin1,vafmin1-vafmin1(1)-vaf0,'linewidth',2);
	plot(tmax1,vafmax1-vafmax1(1)-vaf0,'linewidth',2);
	legend('<1% neg. gradient','>99% pos. gradient',...
		'min grad element','max grad element','location','northwest');
	xlabel('years');
	ylabel('VAF anomaly (compared to n=3) (m^3)'); % m^3


end % }}}
