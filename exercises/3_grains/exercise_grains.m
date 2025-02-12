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







%% 2) make maps of aspect ratio and long axis direction for each phase
% useful commands: smooth, principalComponents, isBoundary, isIndexed








%% 3) which phase shows the strongest grain shape anisotropy?
% useful commands: principalComponents, histogram, eig, paror






%% 4) which crystal directions align parallel to the
%  hornblende grain long axes





%% 5) a) make a colorcoding which highlights
% the alignment of [001], (0,1,0) and [0,0,1] of Hornblende







%%  6) display grains which are more/less than 25*degree with [001] away from
% the [001] maximum







