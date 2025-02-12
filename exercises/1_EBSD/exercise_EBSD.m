% 
% 
% 
%    /$$$$$$$$ /$$$$$$$   /$$$$$$  /$$$$$$$ 
%   | $$_____/| $$__  $$ /$$__  $$| $$__  $$
%   | $$      | $$  \ $$| $$  \__/| $$  \ $$
%   | $$$$$   | $$$$$$$ |  $$$$$$ | $$  | $$
%   | $$__/   | $$__  $$ \____  $$| $$  | $$
%   | $$      | $$  \ $$ /$$  \ $$| $$  | $$
%   | $$$$$$$$| $$$$$$$/|  $$$$$$/| $$$$$$$/
%   |________/|_______/  \______/ |_______/ 
% 
% 
% 
%       /$$$$$$   /$$$$$$   /$$$$$$  /$$   /$$
%      /$$__  $$ /$$$_  $$ /$$__  $$| $$  | $$
%      |__/  \ $$| $$$$\ $$|__/  \ $$| $$  | $$
%       /$$$$$$/| $$ $$ $$  /$$$$$$/| $$$$$$$$
%      /$$____/ | $$\ $$$$ /$$____/ |_____  $$
%      | $$      | $$ \ $$$| $$            | $$
%      | $$$$$$$$|  $$$$$$/| $$$$$$$$      | $$
%      |________/ \______/ |________/      |__/
% 
%
%
%                                    


%% import some data

ebsd = EBSD.load('Rialb_quarzite.h5oina','interface','h5oina','fullDataset')
% note: here we specify 'fullDataset' since the h5oina file for the
% exercise is quite stripped down and some properties were removed to save
% space

% correct the data for inconsistent coordinate systems
% more on that in the next days
rot = rotation.byAxisAngle(yvector,pi);
ebsd  =rotate(ebsd,rot,"keepXY")

% set the plotting convention
pC = plottingConvention(-zvector,xvector);
setMTEXpref('xyzPlotting',pC);
%% 1) inspect the ebsd object



%% 2) plot the image quality with grey colormap and colorbar
% useful commands: plot, mtexColorMap, mtexColorbar



%% 3) select points with a certain property 
% a) plot the mad and plot only those points with 
% a mad < 1 and indexed phase



%% 4) remove and fill (partially)isolated non-indexed points
% useful commands: erode, fill



%% 5) make an ipf colormap
% make an ipf color mapping  with respect to
% x, y,z and a direction of your choice
% useful commands: ipfHSVKey

% ebsd('q').CS=ebsd('q').CS.Laue



%% 6) make a c-axis mapping
% make a c-axis colorcoding and plot a c-axis mapping
% useful commands: makeCaxisDirections, HSVDirectionKey



%% 7) colorize all points 
% as a function of angular distance from the screen normal 
% useful commands: angle



%% 8) quantify the percentage of pixels with c-axes
% close to the screen. Plot maps of those pixels (and those which are not)



%% 8) plot the orientation gradient with respect 
%  to x, y and the average magnitude (in angle/µm)


%%
function cdir = makeCaxisDirections(ebsd)
h = Miller(0,0,0,1,ebsd('q').CS);
cdir = ebsd('q').orientations.*h(1);
% cdir(cdir.z<0) = -cdir(cdir.z<0);
cdir.antipodal=1;
end


