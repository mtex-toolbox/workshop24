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


%% 1c) plot a sharp map of mis2mean
% useful commands: axisAngleColorKey


%% 1d) de-noise the data (optional)
% useful commands: halfQuadraticFilter, smooth

%% 2) inspect grain boundary misorientations





%% 3) find the twinning relation
% useful commands: calcDensity, max, round(...'max')



%% 4) merge the twins
% useful commands: merge, isTwinning











% update the grainIds to the parentIds
ebsd_m('indexed').grainId = parentId(grains.id2ind(ebsd('indexed').grainId))


