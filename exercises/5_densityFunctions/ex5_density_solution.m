%%                        MTEX Workshop 2024
%
%%               Exercises 5 - Density Functions
%
%% 0. Data Import
% Import the data set from exercise 1 and generate an IPF-Z colored map.

ebsd = EBSD.load('Rialb_quarzite.h5oina','interface','h5oina','fullDataset');

% correct the data for inconsistent coordinate systems
rot = rotation.byAxisAngle(yvector,pi);
ebsd = rotate(ebsd,rot,"keepXY");

% set the plotting convention
plottingConvention.default.east = xvector;
plottingConvention.default.outOfScreen = -zvector;

plot(ebsd,ebsd.orientations)

%% 1. ODF
% 
%% a) inverse pole figures
% Extract all Quartz orientations from the ebsd map and visualize them in a
% pole figures corresponding to the planes (0001), (10-10), (10-11)

cs = ebsd.CS;
ori = ebsd.orientations;

h = Miller({0,0,0,1},{1,0,-1,0},{1,0,-1,1},cs)
plotPDF(ori,h)

%% b)
% Visualize the orientation in sigma sections corresponding to the angles 
% sigma = 0:15:105 degree.
% * Try to interpret the orientation data
% * Use the option 'contourf' to generate a smooth plot.
% * Determine the modal orientation by clicking on it with a mouse and
% define it as an orientation ori0.

figure(1)
plotSection(ebsd.orientations,'sigma',(0:15:105)*degree,'contourf')

%%

ori0 = orientation.byEuler(0,30*degree,70*degree,cs)

%% c)
% Compute an ODF from the orientations and visualize it as sigma sections
% * Compare it to the section plot generated in b)
%
% useful commands: orientation/calcDensity

odf = calcDensity(ori)

figure(2)
plotSection(odf,'sigma',(0:15:105)*degree)

%% d)
% Compute ODFs with the option halfwidth set to different values between
% 2.5*degree and 20*degree.
% * What do you observe?
% * What do you think is the best value?

odf = calcDensity(ori,'halfwidth',5*degree)

figure(2)
plotSection(odf,'sigma',(0:15:105)*degree)

%% e)
% Keep going with the odf ODF with best halfwidth value.
% * Evaluate the ODF at ori0 as found in 1b)
% * Determine maximum value of ODF and the corresponding orientation ori1.
% * Determine the misorientation angle to the orientation ori0 found in 1b)
%
% useful commands: SO3Fun/eval, SO3Fun/max

odf.eval(ori0)

%%

[value, ori1] = max(odf)

%%

angle(ori0,ori1) ./ degree

%% f) 
% Compute the volume percentage of EBSD pixels have an orientation ori1
% with a threshold of 10 degree
% * using the the EBSD orientations
% * using the ODF
% * how does the ODF volume changes with the halfwidth
% * compare with the volume of an uniform ODF
%
% useful commands: SO3Fun/volume

v0 = volume(ori,ori1,10*degree) * 100
v1 = volume(odf,ori1,10*degree) * 100

%% e)
% Determine the second local maximum in the ODF 
%

[value,ori1] = max(odf,'numLocal',2)

plotSection(odf,'sigma',(0:10:110)*degree)
annotate(ori1)


%% f)
% Use the command SO3Fun/calcComponents to find the two major texture
% components. Compare the volumes determined by this function with the
% previously found volumes.
% 

[ori2,vol] = calcComponents(odf)

%%
plotSection(odf,'sigma',(0:10:110)*degree)
annotate(ori1(1:2))

%% 2.  Pole Figures

%% a) 
% Use the ODF from exercise 1 to compute the c-axis pole figure of the
% texture.
% * plot it and compare it to the pole figure plotted directly from the
% orientations.
%
% useful commands: SO3Fun/calcPDF

pdf = calcPDF(odf,Miller(0,0,0,1,cs))

plot(pdf)

%% b)
% * Determine the two major maxima in the c-axis pole figure.
% * Compare with the values we would guess from the two modal orientations
% of the texture.


%% c)
% Determine the volume of all grains that are aligned with the c-axis in
% direction of these two maxima within a tolerance of 10 degree.
% * perform the computation with the EBSD data
% * perform the computation with the ODF
%
% useful command: SO3Fun/volume


%% d) 
% Rotate the ODF such that the major c-axis points into z-direction.


%% e)
% Compute a random sample of this c-axis pole density function and check
% how close they are to the original sampling.



%% 3. XRD Data
%
%% a)
% Import the Dubna Neutron test data and inspect them visually.

mtexdata dubna

%% b)
% In each pole figure there are two theta angles for which the intensities
% do not fit very well to the other theta angles. Identify these circles
% and remove them from the data.
%
% useful commands: pf.r, pf.r.theta


%% c)
% Reconstruct an ODF from the neutron data and visualize compare its pole
% density functions with the measured data.
%
% useful commands: ODF/plotPDF


%% d)
% Plot the c-axis pole figure from the reconstructed ODF. Compute and mark
% the preferred direction of the c-axis.


%% e)
% Determine volume percentage of c-axes orientations within a 10 degree
% ball round the preferred c-axis direction. Compare the percentage with an
% untextured ODF.

%% f)
% Visualize the ODF in sigma sections and compare these with the c-axis
% pole figure. 


%% g)
% Find two preferred orientations within the ODF.
%
% useful commands: ODF/max


%% h)
% Try to approximate the reconstructed ODF by an model ODF consisting of 3
% components only. 
%
% useful commands: unimodalODF, uniformODF, fibreODF, binghamODF, 


%% i)
% Determine the error between the ODF computed from the neutron data and
% your approximation in h).
%
% useful commands: ODF/calcError
