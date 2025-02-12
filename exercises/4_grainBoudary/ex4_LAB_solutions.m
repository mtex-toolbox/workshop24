%**************************************************************************
% %**************************************************************************
% 
%
%
%            /$$        /$$$$$$  /$$$$$$$ 
%           | $$       /$$__  $$| $$__  $$
%           | $$      | $$  \ $$| $$  \ $$
%           | $$      | $$$$$$$$| $$$$$$$ 
%           | $$      | $$__  $$| $$__  $$
%           | $$      | $$  | $$| $$  \ $$
%           | $$$$$$$$| $$  | $$| $$$$$$$/
%           |________/|__/  |__/|_______/ 
%                               
%        $$$$$$\   $$$$$$\   $$$$$$\  $$\   $$\
%       $$  __$$\ $$$ __$$\ $$  __$$\ $$ |  $$ |
%       \__/  $$ |$$$$\ $$ |\__/  $$ |$$ |  $$ |
%        $$$$$$  |$$\$$\$$ | $$$$$$  |$$$$$$$$ |
%       $$  ____/ $$ \$$$$ |$$  ____/ \_____$$ |
%       $$ |      $$ |\$$$ |$$ |            $$ |
%       $$$$$$$$\ \$$$$$$  /$$$$$$$$\       $$ |
%       \________| \______/ \________|      \__|
%                           
%                      
%**************************************************************************
%**************************************************************************
    

% Exercise: grain boundaries and misorientations
% we will load a crop from a larger map, previously saved in mtex

load('QvR_8_smoo.mat')
load('Tnl_2_smoo.mat')

%% 1) Compute grains 
% and determine the total grain boundary length. How much does it change
% when smoothign grains?
% useful commands: calcGrains, smooth

[grains ebsd.grainId ebsd.mis2mean] = ebsd.calcGrains('angle',[10 1]*degree);
l1 = sum(grains.boundary('q','q').segLength)

grains = smooth(grains,9,'moveTriplePoints')
l2 = sum(grains.boundary('q','q').segLength)

% change
(l1-l2)/l1


%% 2) Determine low angle (LA) boundaries
% if we computed grains with a second argument (e.g. 1 degree), we have
% LAB down to 1 degree, however they go up to... high angles
% useful commands: innerBoundary

gB = grains.boundary('q','q')
iB = grains('q').innerBoundary

plot(ebsd('q'),ebsd('q').orientations)
hold on 
plot(gB)
hold on
plot(iB,'linecolor','b') %innerBoundary 
hold off

% max angle of inner boundaries
maxiBa = max(grains('q').innerBoundary.misorientation.angle)/degree

plot(ebsd('q'),ebsd('q').orientations,'faceAlpha',0.3)
hold on 
plot(gB)
hold on
plot(iB,iB.misorientation.angle/degree) %innerBoundary 
hold off


% condition for misorientation angles of inner boundaries smaller ...
cond = iB.misorientation.angle < 9*degree;

% those are the low angle boundaries now
LAB = iB(cond)

plot(ebsd('q'),ebsd('q').orientations,'FaceAlpha',0.3)
hold on
plot(gB, 'linecolor','k')
plot(LAB,'edgeAlpha',LAB.misorientation.angle/degree / 9, 'linecolor','k')
hold off


% this is the misorientation of the low angle boundaries
mori = LAB.misorientation;

%% --------------------------------------------------------------------------
% a) Determine the distribution of low angle boundaries (1-9 degree).
% useful commands: plotAngleDistribution plotAxisDistribution
plotAngleDistribution(mori)

plotAxisDistribution(mori)

%% --------------------------------------------------------------------------
% b) Colorcode misorientation axes in crystal coordiantes and plot them in
% a map.
% useful commands: axis, HSVDirectionKey

% misorientation axes in crystal coordinates
axC =  mori.axis;
% colorcoding
cs = ebsd('q').CS;
cKC = HSVDirectionKey(cs)


figure
plot(ebsd('q'),ebsd('q').orientations,'FaceAlpha',0.3)
hold on
plot(LAB,cKC.direction2color(axC),'linewidth',2,'smooth')
hold off
figure
plot(cKC)



%% --------------------------------------------------------------------------
% c) Colorcode misorientation axes in specimen coordiantes and plot them in
%    a map.
% useful commands: axis, HSVDirectionKey


% extract orientations from both sides of the LAB
o = ebsd('id',LAB.ebsdId).orientations

% compute the misorientation axis in specimen coordinates
axS = axis(o(:,1),o(:,2),'antipodal')
% colorcoding
cKS= HSVDirectionKey(specimenSymmetry('-1'))


figure
plot(ebsd('q'),ebsd('q').orientations,'FaceAlpha',0.3)
hold on
plot(LAB,cKS.direction2color(axS),'linewidth',2)
hold off
figure
plot(cKS)



%% --------------------------------------------------------------------------
% d) Plot the misorientation axes in crystal coordiantes in a point plot
%    and colorcode them according to the corresponding direction in
%    specimen coordinates.
%    AND
%    plot the misorientation axes in specimen coordinates and color code
%    with respect to the misorientatin axis in crytsal coordinates.
%
plot(axC,cKS.direction2color(axS),'MarkerSize',4,'fundamentalRegion')
nextAxis
plot(cKS)
nextAxis
plot(axS,cKC.direction2color(axC),'antipodal')
nextAxis
plot(cKC)


%% --------------------------------------------------------------------------
% e) make an axis-angle colorcoding
% useful command: axisAngleColorKey
aAck = axisAngleColorKey(cs,cs,'antipodal')

aAck.maxAngle = 9*degree;

figure
plot(ebsd('q'),ebsd('q').orientations,'FaceAlpha',0.3)
hold on
plot(LAB,aAck.orientation2color(LAB.misorientation),'linewidth',2,'smooth')
hold off
nextAxis
plot(axC,aAck.orientation2color(LAB.misorientation),'fundamentalRegion', ...
     'MarkerSize',4)


%% 4) Select low angle (LA) boundaries
% a) Select all LA boundaries with misorientation axes which are parallel 
%    to [0,0,0,1] (within some threshold).

threshold = 15*degree

refDirC = Miller(0,0,0,1,ebsd('q').CS);

cond = angle(axC,refDirC) < threshold;

LAB(cond);

plot(LAB,'linewidth',2)
hold on
plot(LAB(cond),'LineColor','r','linewidth',2)
hold off
nextAxis

plot(axC(cond).project2FundamentalRegion, ...
     'MarkerFaceColor','r','fundamentalRegion')
hold on
plot(axC(~cond).project2FundamentalRegion, ...
     'MarkerFaceColor','k','fundamentalRegion')
hold on
plot(refDirC.symmetrise('unique'),'MarkerFaceColor','b')
circle(refDirC.symmetrise('unique'),threshold,'linewidth',2,'LineColor','b')
hold off

%% --------------------------------------------------------------------------
% b) Select all LA boundaries with misorientation axes which are parallel to 
%    the z direction. Display the boundaries on a map and on a pole plot.

threshold = 15*degree;
refDirS = zvector;
cond = angle(axS,refDirS) < threshold;
LAB(cond);


plot(LAB,'linewidth',2)
hold on
plot(LAB(cond),'LineColor','r','linewidth',2)
hold off
nextAxis


plot(axS(cond),'MarkerFaceColor','r','antipodal')
hold on
plot(axS(~cond),'MarkerFaceColor','k','antipodal')
hold on
plot(refDirS,'antipodal','MarkerFaceColor','b')
circle(refDirS,threshold,'linewidth',2,'antipodal','LineColor','b')
hold off



%% 5) Can we check if some of these LABs are tilt boundaries
%  A tilt boundary should have the rotation axis within the boundary plane.
%  What is possible with 2D EBSD?

% compute the angle between the grain boundary trace and the misorientation
% axis in specimen coordiantes. Only if those are "parallel" we can be sure
% for this given boudnary that it might be a tilt boundary. Since we lack
% any other 3D information on the grain boudnary, any assumption for tilt
% boundaries would involve some additional info on its inclination.

gbta = angle(axS,LAB.direction);

plot(LAB,gbta/degree,'linewidth',2)
hold on
quiver(LAB(1:5:end),axS(1:5:end),'color','k','autoArrowSize')
hold off




