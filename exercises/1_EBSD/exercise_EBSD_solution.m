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
ebsd

%% 2) plot the image quality with grey colormap and colorbar
% useful commands: plot, mtexColorMap, mtexColorbar

plot(ebsd,ebsd.bc,'showMicronBar','on')
mtexColorMap black2white
mtexColorbar('title','band contrast','Location','southoutside')

%% 3) select points with a certain property 
% a) plot the mad and plot only those points with 
% a mad < 1 and indexed phase

cond1 = ebsd.MAD < 1;
plot(ebsd(cond1),ebsd(cond1).MAD)
nextAxis
cond2 = ebsd.isIndexed;
plot(ebsd(cond1&cond2),ebsd(cond1&cond2).MAD)


%% 4) remove and fill (partially)isolated non-indexed points
% useful commands: erode, fill

ebsd = erode(ebsd,2)
ebsd = fill(ebsd)

%% 5) make an ipf colormap
% make an ipf color mapping  with respect to
% x, y,z and a direction of your choice
% useful commands: ipfHSVKey

% ebsd('q').CS=ebsd('q').CS.Laue

v = [xvector yvector zvector  vector3d.byPolar(pi/2, -30*degree)];
ck = ipfHSVKey(ebsd('q').CS);

for i = 1:length(v)  

    ck.inversePoleFigureDirection=v(i);

    color = ck.orientation2color(ebsd('q').orientations);
    nextAxis
    plot(ebsd('q'),color)
    mtexTitle(['ref. dir. (' char(v(i)) ')'])
   
end

%% 6) make a c-axis mapping
% make a c-axis colorcoding and plot a c-axis mapping
% useful commands: makeCaxisDirections, HSVDirectionKey

cdk=HSVDirectionKey(specimenSymmetry('-1'))

cdir = makeCaxisDirections(ebsd('q'))

color = cdk.direction2color(cdir);

plot(ebsd('q'),color)

%% 7) colorize all points 
% as a function of angular distance from the screen normal 
% useful commands: angle

cdir = makeCaxisDirections(ebsd('q'))

plot(ebsd('q'),angle(zvector,cdir)./degree)


%% 8) quantify the percentage of pixels with c-axes
% close to the screen. Plot maps of those pixels (and those which are not)

cond = angle(zvector,cdir) < 70*degree;
ebsdq = ebsd('q');
plot(ebsdq(cond),'FaceColor',[0.3 0.8 0.2])
hold on
plot(ebsdq(~cond),'FaceColor',[0.9 0.1 0.5])
hold off


%% 8) plot the orientation gradient with respect 
%  to x, y and the average magnitude (in angle/µm)
ebsd  = ebsd.gridify;
gradx = ebsd.gradientX;
grady = ebsd.gradientY;
plot(ebsd,norm(gradx)/degree)
nextAxis
plot(ebsd,norm(grady)/degree)
nextAxis
plot(ebsd,sqrt(norm(gradx).^2 + norm(grady).^2)/degree)

%%
function cdir = makeCaxisDirections(ebsd)
h = Miller(0,0,0,1,ebsd('q').CS);
cdir = ebsd('q').orientations.*h(1);
% cdir(cdir.z<0) = -cdir(cdir.z<0);
cdir.antipodal=1;
end


