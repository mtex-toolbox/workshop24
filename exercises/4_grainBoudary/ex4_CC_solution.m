%**************************************************************************
%**************************************************************************
%                    
%                     $$\           $$\   $$\               
%                     $$ |          \__|  $$ |              
%   $$$$$$$\ $$$$$$\  $$ | $$$$$$$\ $$\ $$$$$$\    $$$$$$\  
%  $$  _____|\____$$\ $$ |$$  _____|$$ |\_$$  _|  $$  __$$\ 
%  $$ /      $$$$$$$ |$$ |$$ /      $$ |  $$ |    $$$$$$$$ |
%  $$ |     $$  __$$ |$$ |$$ |      $$ |  $$ |$$\ $$   ____|
%  \$$$$$$$\\$$$$$$$ |$$ |\$$$$$$$\ $$ |  \$$$$  |\$$$$$$$\ 
%   \_______|\_______|\__| \_______|\__|   \____/  \_______|
%                                                          
%                                                                                                                                  
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
%**************************************************************************
%**************************************************************************
    

% we will load a crop from a larger map, previously saved in mtex

% ebsd = loadEBSD_h5oina('Juf19nN_ss1.h5oina')
%   
% 
% % correct the data for inconsistent coordinate systems
% rot = rotation.byAxisAngle(yvector,pi);
% ebsd  =rotate(ebsd,rot,"keepXY")

% or 
 ebsd = EBSD.load('Juf19nN_ss1.ctf')
% no correction needed with this one as it has been exported from mtex


%% 1a) Compute grains, get rid of small, non-indexed grains
% and make smooth outlines
% useful commands: calcGrains, smooth

[grains ebsd.grainId ebsd.mis2mean] = ebsd.calcGrains('angle',[10 1]*degree);
ebsd(grains(grains.grainSize<6 & ~grains.isIndexed))=[];
[grains ebsd.grainId ebsd.mis2mean] = ebsd.calcGrains('angle',[10 1]*degree);


grains = smooth(grains,6,'moveTriplePoints')

plot(grains)

%% 1b) plot an ipf map with respct to x,y,z...
% useful commands: ipfHSVKey
cs = ebsd('c').CS
ck =ipfHSVKey(cs)
ck.inversePoleFigureDirection=xvector
plot(ebsd('c'),ck.orientation2color(ebsd('c').orientations))
hold on
plot(grains.boundary)
hold off


%% 1c) plot a sharp map of mis2mean
% useful commands: axisAngleColorKey
ck3 = axisAngleColorKey(cs,cs);
ck3.maxAngle = 3*degree;
plot(ebsd('c'),ck3.orientation2color(ebsd('c').mis2mean))


%% 1d) de-noise the data
% useful commands: halfQuadraticFilter, smooth

%set up a filter
hq = halfQuadraticFilter
hq.alpha = 0.75;
hq.threshold = 10*degree;

%also; let's fill the single/double nonIndexed parts
ebsd = smooth(ebsd,hq,'fill')

% compute only mis2mean again (why is that not so smart?)
ebsd.mis2mean=calcGROD(ebsd,grains)
nextAxis
plot(ebsd('c'),ck3.orientation2color(ebsd('c').mis2mean))

% compute grains again
[grains ebsd.grainId ebsd.mis2mean] = ebsd.calcGrains('angle',[10 1]*degree);
nextAxis
plot(ebsd('c'),ck3.orientation2color(ebsd('c').mis2mean))


%% 2) inspect grain boundary misorientations
gbcc = grains.boundary('c','c')

mori = gbcc.misorientation;

figure
plot(ebsd('c'),ebsd('c').orientations,'FaceAlpha',0.3)
hold on
plot(gbcc,mori.angle/degree,'linewidth',2,'smooth')
hold off
mtexColorbar
%% inspect misorientation axes

axC = mori.axis;

% colorcoding
cs = ebsd('c').CS;
cKC = HSVDirectionKey(cs)


plot(ebsd('c'),ck.orientation2color(ebsd('c').orientations),'FaceAlpha',0.3)
hold on
plot(gbcc,cKC.direction2color(axC),'linewidth',1,'smooth')
hold off
%%
figure
plot(cKC)


%% find the twinning relation
% useful commands: calcDensity, max, round(...'max')


% discard small angles
cond = mori.angle>10*degree;

mdf = calcDensity(mori(cond),'halfwidth',5*degree)

plot(mdf,'contourf','axisAngle')
%
[~,mmori]=max(mdf)
hold on
plot(mmori,'MarkerFaceColor','r','axisAngle')
hold off


round(mmori.axis('max'))
round(mmori.axis)

mmori.angle('max')/degree
mmori.angle/degree

%% merge the twins
% useful commands: merge

gbc = grains.boundary;
ind = gbc.isTwinning(mmori,10*degree);

plot(ebsd('c'),ebsd('c').orientations,'FaceAlpha',0.3)
hold on
plot(gbc(ind))
hold on

[grains_m,parentId] = merge(grains,gbc(ind));
plot(grains_m.boundary('c','c'),'linewidth',2)
hold off


ebsd_m = ebsd;

% update the grainIds to the parentIds
ebsd_m('indexed').grainId = parentId(grains.id2ind(ebsd('indexed').grainId))


% find and plot the largest grain
gid = grains_m(grains_m.grainSize==max(grains_m.grainSize)).id
    
plot(grains_m(gid))
hold on
plot(ebsd_m(grains_m(gid)))
hold off




