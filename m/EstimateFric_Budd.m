function fric=EstimateFric_Budd(md) % {{{

   disp('Estimating friction coefficient - Budd')
   % initial guess from driving stress
   min_velocity=0.1; % m/yr
   min_pressure=1.;% Pa
   min_c=1;

   asurf=averaging(md,md.geometry.surface,20); % maybe executing 20 L2 projection is ok
   [sx,sy,s]=slope(md,asurf); % slope 's' comes on elements
   sslope=averaging(md,s,1); % average the slope once on the vertices, because 's' comes on elements, we need this data on vertices

   vel=md.inversion.vel_obs;
   vel=max(vel,min_velocity); % setting minimum velocity value
   vel=vel/md.constants.yts; % S.I.

   Neff=(md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base)*md.constants.g;
   Neff(find(Neff<=0))=min_pressure; % setting minimum positve pressure

   driving_stress=md.materials.rho_ice*md.constants.g*md.geometry.thickness.*(sslope);
   c=sqrt(driving_stress./(Neff.*vel));
   c=max(c,min_c);

   fric=c; % initial guess

end%}}}
