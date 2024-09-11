function [strainrate] = strainrate_SSA(md,vx,vy,average)
%STRAINRATE_SSA returns the strain rate matrix and effective strain rate using SSA
% As per SSA, the XZ and YZ components are assumed to be zero. 
% As per the incompressibility assumption, we set ZZ = - XX - YY
%
%INPUT:
% 	md		the model you are using
% 	vx,vy		the velocity field components from which to calculate the strainrate, i.e.
%				md.initialization.vx/vy				-the observed velocities
% 				md.results.StressbalanceSolution.Vx/Vy		-the modeled stress balance velocities
%  average  either 1 (default, interpolate back onto each node) or 0 (return as element wise)
%OUTPUT:
% 	strainrate		structure with fields
% 	  xx,xy,yy,zz		each contains the respective component of the 3D SSA strain rate tensor for each node 
% 	  eff			the effective strain rate for each node
%
%USE:
% 	[strainrate] = strainrate_SSA(md,md.initialization.vx,md.initialization.vy); % return node-wise
% 	[strainrate] = strainrate_SSA(md,md.initialization.vx,md.initialization.vy,0); % return element-wise
%
%Benjamin Getraer, main code adapted from Mathieu Morlighem
%Written: 11/2/2022
%Last edited: 11/8/2022

if nargin ==3
	average = 1;
end

%Compute strain rate (adapted from Mathieu)
%*********************************
% extract numbering index for each element
index = md.mesh.elements; 
% extract the nodal coefficients for each element's shape function
[alpha beta]=GetNodalFunctionsCoeff(index,md.mesh.x,md.mesh.y);
% extract the x and y components of the velocity at every node (3 nodes per element)
vxlist = vx(index)/md.constants.yts;
vylist = vy(index)/md.constants.yts;
% sum the component of the gradient of the velocity field over the nodes of each element 
%		alpha are the x derivatives, beta are the y derivatives
ux=sum(vxlist.*alpha, 2);	% ux = eps_xx
uy=sum(vxlist.*beta , 2);	% eps_xy = 1/2(uy+vx) = eps_yx
vx=sum(vylist.*alpha, 2);	% eps_xy = 1/2(uy+vx) = eps_yx	
vy=sum(vylist.*beta , 2);	% vy = eps_yy
if average
	% build the strain rate components and interpolate from the elements back onto each node
	strainrate.xx = averaging(md,ux,0)*md.constants.yts; %strain rate in 1/a instead of 1/s
	strainrate.xy = averaging(md,1/2*(uy+vx),0)*md.constants.yts;
	strainrate.yy = averaging(md,vy,0)*md.constants.yts; 
	strainrate.zz = averaging(md,-ux-vy,0)*md.constants.yts;
else
	% build the strain rate components element-wise
	strainrate.xx = ux*md.constants.yts; %strain rate in 1/a instead of 1/s
	strainrate.xy = 1/2*(uy+vx)*md.constants.yts;
	strainrate.yy = vy*md.constants.yts; 
	strainrate.zz = (-ux-vy)*md.constants.yts;
end
% solve for the effective strain rate (the II invariant)
strainrate.eff = sqrt(strainrate.xx.^2 + strainrate.yy.^2 + strainrate.xy.^2 + strainrate.xx.*strainrate.yy);
