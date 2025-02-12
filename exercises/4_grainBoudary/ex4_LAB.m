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
% load('Tnl_2_smoo.mat')

%% 1) Compute grains 
% and determine the total grain boundary length. How much does it change
% when smoothign grains?
% useful commands: calcGrains, smooth




%% 2) Determine low angle (LA) boundaries
% if we computed grains with a second argument (e.g. 1 degree), we have
% LAB down to 1 degree, however they go up to... high angles
% useful commands: innerBoundary



%% --------------------------------------------------------------------------
% a) Determine the distribution of low angle boundaries (1-9 degree).
% useful commands: plotAngleDistribution plotAxisDistribution



%% --------------------------------------------------------------------------
% b) Colorcode misorientation axes in crystal coordiantes and plot them in
% a map.
% useful commands: axis, HSVDirectionKey




%% --------------------------------------------------------------------------
% c) Colorcode misorientation axes in specimen coordiantes and plot them in
%    a map.
% useful commands: axis, HSVDirectionKey




%% --------------------------------------------------------------------------
% d) Plot the misorientation axes in crystal coordiantes in a point plot
%    and colorcode them according to the corresponding direction in
%    specimen coordinates.
%    AND
%    plot the misorientation axes in specimen coordinates and color code
%    with respect to the misorientatin axis in crytsal coordinates.



%% --------------------------------------------------------------------------
% e) make an axis-angle colorcoding
% useful command: axisAngleColorKey(cs,cs,'antipodal')


%% 4) Select low angle (LA) boundaries
% a) Select all LA boundaries with misorientation axes which are parallel 
%    to [0,0,0,1] (within some threshold).



%% --------------------------------------------------------------------------
% b) Select all LA boundaries with misorientation axes which are parallel to 
%    the z direction. Display the boundaries on a map and on a pole plot.




%% 5) Can we check if some of these LABs are tilt boundaries
%  A tilt boundary should have the rotation axis within the boundary plane.
%  What is possible with 2D EBSD?

% compute the angle between the grain boundary trace and the misorientation
% axis in specimen coordiantes. Only if those are "parallel" we can be sure
% for this given boudnary that it might be a tilt boundary. Since we lack
% any other 3D information on the grain boudnary, any assumption for tilt
% boundaries would involve some additional info on its inclination.



