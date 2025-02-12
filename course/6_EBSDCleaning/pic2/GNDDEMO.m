%% Dislocation Density Estimation
% 
% This example sheet describes how to estimate dislocation densities in
% MTEX following the reference paper
%
% Pantleon, Resolving the geometrically necessary dislocation content by
% conventional electron backscattering diffraction, Scripta Materialia,
% 2008
%


%% data import and grain reconstruction
% Lets start by importing orientation data from 2 percent uniaxial deformed
% steel DC06. 

% set up the plotting convention
plotx2east

% import the EBSD data
%ebsd = loadEBSD('DC06_2uniax.ang');
%ebsd = loadEBSD('DC06_2biax.ang');
ebsd = loadEBSD('DC06_undeformed.ang');

ebsd = rotate(ebsd,90*degree);

% reconstruct grains
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',5*degree);

% remove small grains
ebsd(grains(grains.grainSize<=5)) = [];

% redo grain reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',2.5*degree);

% smooth grain boundaries
grains = smooth(grains,5);

% denoise orientation data
F = halfQuadraticFilter;
F.threshold = 1.5*degree;
F.eps = 1e-2;
F.alpha = 0.005;
ebsd = smooth(ebsd('indexed'),F,'fill',grains);

% consider only the Fe (alpha) phase 
ebsd = ebsd('indexed').gridify

%% plot data

ipfKey = ipfHSVKey(ebsd)
ipfKey.inversePoleFigureDirection = yvector
plot(ebsd,ipfKey.orientation2color(ebsd.orientations),'micronBar','off','figSize','medium')

hold on
plot(grains.boundary,'linewidth',2)
hold off

%saveFigure('uni_ebsd.png')

%% The incomplete curvature tensor
% For gridified orientation maps the command <EBSDSquare.curvature.html
% curvature> computes the incomplete curvature tensor for all positions of
% the map except for the grain boundaries, where the value is set to NaN.
% The unit of the curvature tensor is one over the unit of the spatial data
% in the EBSD map. In the present case this is 1/um as stored inside the
% tensor field.

% compute the curvature tensor
kappa = ebsd.curvature

% one can index the curvature tensors in the same way as the EBSD data.
% E.g. the curvature in pixel (2,3) is
kappa(2,3)

%% The components of the curvature tensor
% As expected the curvature tensor is NaN in the third column as this
% column corresponds to the partial derivative in z-direction which is
% usually unknown for 2d EBSD maps. 
%
% We can access the different components of the curvature tensor with

kappa12 = kappa{1,2};

size(kappa12)

%%
% which results in a variable of the same size as our EBSD map. This allows
% us to visualize the different components of the curvature tensor

mtexFig = newMtexFigure('nrows',3,'ncols',3,'figSize','huge');

plot(ebsd,kappa{1,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{1,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{1,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{2,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{2,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{2,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{3,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{3,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{3,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

% unify the color rage  - you may also use setColoRange equal
setColorRange([-0.005,0.005])
drawNow(gcm,'figSize','large')

%saveFigure('uni_curvatureIncomplete.png')

%% The incomplete dislocation density tensor
% According to Kröner the curvature tensor is directly related to the
% dislocation density tensor. 

alpha = kappa.dislocationDensity

%%
% which has the same unit as the curvature tensor and is incomplete as well
% as we can see when looking at a particular one.

alpha(2,3)

%% Fitting Dislocations to the incomplete dislocation density tensor
% The central idea of Pantleon is that the dislocation density tensor is
% build up by dislocations with different densities such that the total
% energy is minimum. Depending on the attomic lattice different
% dislocattion systems has to be considered. In present case of a body
% centered cubic (bcc) material 48 edge dislocations and 4 screw
% dislocations have to be considered. Those principle dislocations are
% defined in MTEX either by their Burgers and line vectors or by

dS = dislocationSystem.bcc(ebsd.CS)

%%
% Here the norm of the Burgers vectors is important

% size of the unit cell
a = norm(ebsd.CS.aAxis)

% for edge dislocations in bcc the norm of the burgers vector is sqrt(3)/2 * a
[norm(dS(1).b), sqrt(3)/2 * a]

% for screw dislocations it is sqrt(3) * a
% Wolfgang: is this correct or is the norm of the Burgers vector for screw
% dislocations the same as for edge dislocations
[norm(dS(end).b), sqrt(3) * a]

%%
% Note that the energy of each dislocation system can be stored in the
% property |u|

% energy of the edge dislocations
dS(dS.isEdge).u = 1;

% energy of the screw dislocations
dS(dS.isScrew).u = 1 - 0.3;

% Wolfgang: How to do this best? I found different formulas:
%
% E = 1 - poisson ratio
% E = c * G * |b|^2,  - G - Schubmodul / Shear Modulus Energy per (unit length)^2

%%
% A single dislocation causes a deformation that can be represented by a
% the rank one tensor

dS(1).tensor

%%
% Note that the unit of this tensors is the same as the unit used for
% describing the length of the unit cell, which is in most cases Angstrom
% (au). Furthremore, we observe that the tensor is given with respect to
% the crystal reference frame while the dislocation densitiy tensors are
% given with respect to the specimen reference frame. Hence, to make them
% compatible we have to rotate the dislocation tensors into the specimen
% reference frame as well. This is done by

dSRot = ebsd.orientations * dS


%% Soving the fitting problem
% Now we are ready for fitting the dislocation tensors to the dislocation
% densitiy tensor in each pixel of the map. This is done by the command
% <curvatureTensor.fitDislocationSystems.html fitDislocationSystems>.

[rho,factor] = fitDislocationSystems(kappa,dSRot);

%%
% As result we obtain a matrix of densities |rho| such that the product
% with the dislocation systems yields the incomplete dislocation density
% tensors derived from the curvature, i.e.,

% the restored dislocation density tensors 
alpha = sum(dSRot.tensor .* rho,2);

% we have to set the unit manualy since it is not stored in rho
alpha.opt.unit = '1/um';

% the restored dislocation density tensor for pixel 2
alpha(2)

% the dislocation density dervied from the curvature in pixel 2
kappa(2).dislocationDensity

%%
% we may also restore the complete curvature tensor with

kappa = alpha.curvature 

%%
% and plot it as we did before

mtexFig = newMtexFigure('nrows',3,'ncols',3,'figSize','huge');

plot(ebsd,kappa{1,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{1,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{1,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{2,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{2,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{2,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{3,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{3,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

nextAxis
plot(ebsd,kappa{3,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',2); hold off

drawNow(gcm,'figSize','large')
setColorRange([-0.005,0.005])

%saveFigure('uni_curvatureComplete.png')

%% The total dislocation energy 
% The unit of the densities |h| in our example is 1/um * 1/au where 1/um
% comes from the unit of the curvature tensor an 1/au from the unit of the
% Burgers vector. In order to transform |h| to SI units, i.e., 1/m^2 we
% have to multiply it with 10^16. This is exactly the values returned as
% the second output |factor| by the function
% <curvatureTensor.fitDislocationSystems.html fitDislocationSystems>.
  
factor

%% 
% Multiplying the densities |rho| with this factor and the individual
% energies of the the dislocation systems we end up with the total
% dislocation energy. Lets plot this at a logarithmic scale

close all
plot(ebsd,factor*sum(abs(rho .* dSRot.u),2),'micronbar','off')
%plot(ebsd,log(factor*sum(abs(rho .* dSRot.u),2))/log(10),'micronbar','off','figsize','medium')
mtexColorMap('parula')
mtexColorbar

set(gca,'ColorScale','log');
set(gca,'CLim',[1e11 5e14]);
%set(gca,'CLim',[5e10 8e13]);
%set(gca,'CLim',[10 14]);
%set(gca,'CLim',[12 14]);

hold on
plot(grains.boundary,'linewidth',2)
hold off

%saveFigure('gnd.pdf')