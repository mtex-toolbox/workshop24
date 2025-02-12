%==========================================================================
%==========================================================================
%
%                                 $$\
%                                 \__|
%     $$$$$$\   $$$$$$\  $$$$$$\  $$\ $$$$$$$\   $$$$$$$\
%    $$  __$$\ $$  __$$\ \____$$\ $$ |$$  __$$\ $$  _____|
%    $$ /  $$ |$$ |  \__|$$$$$$$ |$$ |$$ |  $$ |\$$$$$$\
%    $$ |  $$ |$$ |     $$  __$$ |$$ |$$ |  $$ | \____$$\
%    \$$$$$$$ |$$ |     \$$$$$$$ |$$ |$$ |  $$ |$$$$$$$  |
%     \____$$ |\__|      \_______|\__|\__|  \__|\_______/
%    $$\   $$ |
%    \$$$$$$  |
%     \______/
%           $$$$$$\   $$$$$$\   $$$$$$\  $$\   $$\
%          $$  __$$\ $$$ __$$\ $$  __$$\ $$ |  $$ |
%          \__/  $$ |$$$$\ $$ |\__/  $$ |$$ |  $$ |
%          $$$$$$  |$$\$$\$$ | $$$$$$  |$$$$$$$$ |
%          $$  ____/ $$ \$$$$ |$$  ____/ \_____$$ |
%          $$ |      $$ |\$$$ |$$ |            $$ |
%          $$$$$$$$\ \$$$$$$  /$$$$$$$$\       $$ |
%          \________| \______/ \________|      \__|
%
%
%
%
%==========================================================================
%==========================================================================

% create an EBSD variable containing the data
% fname = 'SK20_3_3.ctf';
% ebsd = EBSD.load([inpath fname]);
ebsd = loadEBSD_h5oina('Sk20_3_LAM2_cropped.h5oina','fullDataset');

% correct the data for inconsistent coordinate systems
rot = rotation.byAxisAngle(yvector,pi);
ebsd  =rotate(ebsd,rot,'keepXY')

% or in case you have trouble with the h5oina
% load('Sk20_3_LAM2_cropped.mat')
%% inspect the data
plot(ebsd,ebsd.bc)
mtexColorMap black2white
hold on
plot(ebsd,'FaceAlpha',0.2)
hold off

%% 0) clean up the data and compute grains
% suggestion:
% a) backup ebsd.prop (optional)
%    ebsd = ebsd.gridify; props = ebsd.prop;
% b) erode isolated indexed pixels and isolated notIndexed pixels (optional)
% c) compute grains and get rid of unindexed area which are i) small
%    large or ii) not holes or other unknown phases as well as very small
%    indexed grains
% useful commands: gridify, erode, calcGrains, grainMean


ebsd = ebsd.gridify;
props = ebsd.prop;

ebsd = erode(ebsd,3);
ebsd = erode(ebsd,1,{'q' 'al' 'an' 'h' 'ti'});
ebsd = erode(ebsd,1,{'q' 'an' 'h' 'ti'});


[grains ebsd.grainId ebsd.mis2mean]= ebsd.calcGrains
plot(grains)

%% compute mean bc for each grain
% assumption: holes will have a low bc while valid phases which just have
% not been indexed might have higher values
gmbc = grainMean(ebsd,ebsd.bc,grains,@(x) quantile(x,0.95));

plot(grains,gmbc)
hold on
plot(grains('n').boundary,'lineColor','w')
hold off
%%
% set up a condition such that
%    (mean bc is small OR grain size is small) AND in both case the grain
%    is not indexed
%    OR it uis simply a very small grain
cond = ((gmbc < 100 | grains.grainSize < 5) & ~grains.isIndexed) | ...
    grains.grainSize < 3;

ebsd(grains(cond))=[];

[grains ebsd.grainId ebsd.mis2mean]= ebsd.calcGrains
plot(grains)

%%
% since we eroded alls sorts of points but maybe we want to keep valid info
% for some, we could write back the props
plot(ebsd,ebsd.prop.Al_Wt_)
ebsd = ebsd.gridify;
ebsd.prop = props;

nextAxis
plot(ebsd,ebsd.prop.Al_Wt_)


%% 2) make maps of aspect ratio and long axis direction for each phase
% useful commands: principalComponents

grains = smooth(grains,5,'moveTriplePoints')

% note: this is different to previous behavior
[a,b] = principalComponents(grains);

asr = norm(a) ./ norm(b);
omega = mod(a.rho,pi);

cond = grains.isIndexed & ~grains.isBoundary;
plot(grains(cond),asr(cond))
mtexColorbar
setColorRange([1 3])

%% aspect ratio maps for each phase
close 
pid  = grains.indexedPhasesId;

for i  = 1:length(pid)
    cond = grains.phaseId==pid(i) & ~grains.isBoundary;
    nextAxis
    plot(grains(cond),asr(cond))
    mtexTitle(grains(cond).CS.mineral)
end
setColorRange([1 3])

%% long axis maps for each phase
figure
for i  = 1:length(pid)
    cond = grains.phaseId==pid(i) & ~grains.isBoundary;
    nextAxis
    plot(grains(cond),omega(cond)/degree)
    mtexTitle(grains(cond).CS.mineral)
end
mtexColorMap hsv
setColorRange([0 180])

%% 3) which phase shows the strongest grain shape anisotropy?
close

% [a,b] = principalComponents(grains);
% asr = norm(a) ./ norm(b);


pid  = grains.indexedPhasesId;

for i  = 1:length(pid)
    cond = grains.phaseId==pid(i) & ~grains.isBoundary & grains.grainSize > 25;

    subplot(1,length(pid),i)

    lax = grains(cond).longAxis;
    asrG = grains(cond).aspectRatio;
    histogram(lax,'weights',asrG,36)

    % histogram(a(cond),'weights',asr(cond),36)

    % compute the eigenvector, maybe weight with the aspect ratio
    [ev,eva]=eig(lax .* asrG);
    % [ev,eva]=eig(a(cond) .* asr(cond));

    title([grains(cond).CS.mineral ' e: ' num2str(round(eva(3)/sum(eva),2))])
end


%% do the paror
figure
omega = 0:180;
for i  = 1:length(pid)
    cond  = grains.phaseId==pid(i) & ~grains.isBoundary  & grains.grainSize > 25;
    cpf   = grains(cond).paror(omega*degree);
    plot(omega,cpf)
    hold on
end
hold off
legend({'q' 'a' 'h' 't'})

%% 4) which crystal directions align parallel to the
%  hornblende grain long axes
[ah,bh] = principalComponents(grains('h'));

% crystal direction parallel to ah
glh = inv(grains('h').meanOrientation) .* ah

% color code this crystal direction on the map
ck = HSVDirectionKey(glh.CS)
plot(grains('h'),ck.direction2color(glh))
nextAxis
plot(glh.project2FundamentalRegion,'fundamentalRegion')


% for rather equiaxes grains, the long axis is poorly defined
asrh = norm(ah) ./ norm(bh);
plot(glh.project2FundamentalRegion,asrh,'fundamentalRegion')
setColorRange([1 max(asrh)])
sf = calcDensity(glh,'weights',asrh)
[~,mpos] = max(sf)
nextAxis
plot(sf,'contourf')
annotate(mpos)
%% 5) a) make a colorcoding which highlights
% the alignment of [001], (0,1,0) and [0,0,1] of Hornblende

h=Miller({0,0,1},{1,0,0},ebsd('h').CS,'uvw')
hk=Miller(0,1,0,ebsd('h').CS)


plotPDF(grains('h').meanOrientation,h)
nextAxis
plotPDF(grains('h').meanOrientation,hk)


% compute the direction of the max in specimen coordinates
cdir = grains('h').meanOrientation.*h(1);
sf = calcDensity(cdir,'weights',grains('h').area)
sf.antipodal=1;
[~,cmax] = max(sf)
plot(sf)
annotate(cmax,'antipodal')
figure
plot(grains('h'),angle(cmax,cdir)/degree)




%%  5) display grains which are more/less than 25*degree with [001] away from
% the [001] maximum
cond = angle(cmax,cdir)>25*degree;
hgrains = grains('h')
plot(hgrains(cond),'FaceColor','Fuchsia','lineColor','k')
hold on
plot(hgrains(~cond),'FaceColor',[0.2 0.5 1],'lineColor','k')
hold off





