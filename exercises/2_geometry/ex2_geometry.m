%%                        MTEX Workshop 2024
%
%%                   Exercise 2 - Crystal Geometry
%
% * The exercises will be discussed during the second exercise session.
% * You might want to go through it before to check what you can do for
% your own.
% * The main purpose of the exercises is to get used to the syntax of MTEX
% and learn some, maybe not so well known commands.
% * Try to make to figure nice. Adjust color, linewidth, ... and ask if you
% want to know about the options.
% * Try to use tab completion as much as possible. It prevents you of
% misspelling commands and options.
% * The exercises contain hints for useful commands. Type

help commandName

% into the command window to get additional information to this command


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

%% 1. Extract orientations
% Extract all orientations of the data set. How can you distinguish between
% Quartz orientations and the orientations corresponding to not indexed
% pixels?

ori = ebsd.orientations

%%

ebsd('notIndexed').id(1:10)

%%

sum(ori.isnan)


%% 2. Visualize the orientations
%
%% 2a 
% Display some of the orientations in three dimensional Euler and axis
% angle space.
% 
% useful commands: orientation/plot

plot(ori)

%%

plot(ori,'axisAngle','points',10000)


%% 2b
% Display some of the orientations in two dimensional sections of the
% orientations space
% 
% useful commands: orientation/plotSection


plotSection(ori,'phi2',90*degree)


%% 3. Pole Figures
%
%% 3a
% Extract and crystal symmetry of quartz and define the quartz crystal
% directions corresponding to the c-axis, the hexagonal prism and the
% positive rhomboedron lattice plane
%
% useful commands: Miller

h = Miller({0,0,0,1},{1,0,-1,0},{1 0 -1 1}  ,ebsd('qu').CS)

hh = h.symmetrise


%% 3b
% Visualize these crystal directions in a Schmidt net.
% * colorize the crystal directions
% * display all symmetrically equivalent crystal directions
% * display the names of the crystal directions in a legend
%
% useful commands: vector3d/scatter, Miller/scatter, hold on,
% 'displayName', legend

for i = 1:3
  plot(h(i).symmetrise('unique'),'upper','DisplayName',char(h(i),'LaTex'),'labeled')
  hold on
end
hold off

leg = legend 
leg.Interpreter = "latex"

%% 3c
% Compute the angles between all positive rhomboedron lattice planes
%
% * give the angles in degree
% * what is the maximum angle?
%
% useful commands: Miller/symmetrise, Miller/angle


angle(h(3),h(3).symmetrise,'noSymmetry') / degree

%%
clf
plot(h(3).symmetrise,'3d')



%% 3d
% Visualize the orientations in the EBSD dataset in pole figures
% corresponding to the three above defined crystal directions
%
% * try the option 'MarkerFaceAlpha' to see where the orientations become
% more dense
% * try the options 'contourf'
% * use your mouse to detect the most prominent c-Axis position and define
% it as a the variable cAxPos of type vector3d
%
% useful commands: orientation/plotPDF vector3d/byPolar

plotPDF(ori,h,'MarkerAlpha',0.01)

%%
clf
plotPDF(ori,h,'pcolor')

%%

cAxPos = vector3d.byPolar(32*degree,-92*degree)

annotate(cAxPos)


%% 4 Inverse Pole Figures
%
%% 4a
% Define the specimen direction x,y,z and visualize the Quartz orientations
% in the corresponding inverse pole figures. 
%
% * Which planes are aligned with the x-axis of the measurement plane
%
% useful commands: orientation/plotIPDF


r = [xvector,yvector,zvector];

plotIPDF(ori, r, 'contourf')



%% 4b
% Visualize the inverse pole figure with respect to cAxPos as found in
% exercise 3d. Explain the result!

plotIPDF(ori, cAxPos,'contourf')


%% 4c
% We want to analyze whether there is a correlation between the hexagonal
% prism plane of the crystals and the x direction.
% 
% * Determine the angles between the prism plane and the x-axis
% * Plot an histogram of these angles in degree with 20 bins
% * Simulate as many random orientations as orientations in the dataset
% * Overlay the histogram with a histogram of the random orientations
%
% useful commands: histogram, orientation.rand, Miller/angle, hold on


hh = inv(ori)*xvector

%%


% angle between the prism plane to specimen x-axis
omega = angle(h(2), hh) / degree;

clf
histogram(omega,180)

%%

oriRand = orientation.rand(length(ori),ori.CS)

omegaRand = angle(h(2), inv(oriRand)*xvector) / degree;

hold on
histogram(omegaRand,180)
hold off


%% 4d
% Next we want to give exact volume portions ....
%
% * Determine the percentage of orientations that have their prism
% plane aligned with the x-axis within a tolerance of 10 degree.
% * compare this value with the percentage of random orientations that have
% this property.
% * What is the quotient between these two values?
%
% useful commands: sum, nnz, length


100 * sum(omega < 10) / sum(~isnan(omega))
100 * sum(omegaRand < 10) / sum(~isnan(omegaRand))


%% 5. Rotate the EBSD Data
%
%% 5a
% Determine a rotation that maps the major c-axis position to z and keeps
% the x direction.
%
% useful commands: rotation.byAxisAngle, rotation.map, rotation.fit


%% 5b
% * Rotate the entire EBSD map by the above rotation.
% * Plot an ipf-z colored map of the rotated data.
% * Compare it to the original ipf-z colored map.


%% 5c
% Compare the pole figures of the rotated EBSD map with the original pole
% figures.



%% 5d
% Compare the inverse pole figures of the rotated EBSD map with the
% original pole figures.



%% 5e
% Alternatively, we could also change the plotting convention such that the
% major c-axis points out of plane and the x-axis is kept pointing east
%
% * define a new plotting convention that does this
% * use this plotting convention for displaying the original pole figures
% * use this plotting convention for displaying the original EBSD map
%



%% 6. Cluster Analysis
% Next we want to determine clusters of dominant orientations.
%
%% 6a
% Determine clusters of similarly oriented pixels.  
%
% useful commands: orientation/calcCluster (takes some time)


%% 6b
% Colorize the orientations according to the cluster they belong to
%
% * in the axis angle plot
% * in the EBSD map
% * in the pole figures
% * in the inverse pole figures



%% 6c Volume portions
% Determine the percentage of orientations that belong to cluster 1 and 2
%


