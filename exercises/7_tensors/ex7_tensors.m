%% MTEX workshop 204
% Exercise of physical properties
% Seismic properties of mantle rock

%% 1. Close all figures and clean workspace

close all
clear all

% Setup setMTEXpref - basiclcc annotations are off
setMTEXpref('FontSize',25);
setMTEXpref('pfAnnotations',@(varargin) 1);

%% 2. Import some data

mtexdata forsterite
ebsd = rotate(ebsd,rotation('axis',-zvector,'angle',90*degree));

%% 2.1 Plot orientation map for forsterite

 


%% 3. Import diopside elastic tensor

cs_Tensor_cpx = crystalSymmetry('121',[9.585  8.776  5.26],...
  [90.0000 105.8600  90.0000]*degree,'x||a*','z||c',...
  'mineral','Diopside');

% the density
rho_cpx = 3.2860;

% the tensor coefficients
Cij_cpx =....
  [[  228.10   78.80   70.20    0.00    7.90    0.00];...
  [    78.80  181.10   61.10    0.00    5.90    0.00];...
  [    70.20   61.10  245.40    0.00   39.70    0.00];...
  [     0.00    0.00    0.00   78.90    0.00    6.40];...
  [     7.90    5.90   39.70    0.00   68.20    0.00];...
  [     0.00    0.00    0.00    6.40    0.00   78.10]];

% define the tensor
C_cpx = stiffnessTensor(Cij_cpx,cs_Tensor_cpx,'density',rho_cpx);

%% 3.1 Plot seismic velocities of diopside single crystal


%% 4. Import olivine elastic tensor

cs_Tensor_ol = crystalSymmetry('222', [4.762 10.225 5.994],...
    'mineral', 'olivine', 'color', 'light red');

% we need the density of this phase
rho_ol = 3.3550;

% and we need the coefficients Cij for olivine in Voigt notation 
% Abramson et al. 1997 Journal of Geophysical Research)

Cij_ol = [[320.5  68.15  71.6     0     0     0];...
       [ 68.15  196.5  76.8     0     0     0];...
       [  71.6   76.8 233.5     0     0     0];...
       [   0      0      0     64     0     0];...
       [   0      0      0      0    77     0];...
       [   0      0      0      0     0  78.7]];

% the elastic stiffness tensor is defined as

C_ol = stiffnessTensor(Cij_ol,cs_Tensor_ol,'density',...
    rho_ol);

%% 4.1 Plot Young Modulus forsterite



%% 5. Import enstatite elastic tensor

cs_Tensor_opx = crystalSymmetry('mmm',[ 18.2457  8.7984  5.1959],...
  [  90.0000  90.0000  90.0000]*degree,'x||a','z||c',...
  'mineral','Enstatite');

% the density
rho_opx = 3.3060;

% the tensor coefficients
Cij_opx =....
  [[  236.90   79.60   63.20    0.00    0.00    0.00];...
  [    79.60  180.50   56.80    0.00    0.00    0.00];...
  [    63.20   56.80  230.40    0.00    0.00    0.00];...
  [     0.00    0.00    0.00   84.30    0.00    0.00];...
  [     0.00    0.00    0.00    0.00   79.40    0.00];...
  [     0.00    0.00    0.00    0.00    0.00   80.10]];

% define the tensor
C_opx = stiffnessTensor(Cij_opx,cs_Tensor_opx,'density',rho_opx);

%% 6. Calculate the ODFs of the individual phases

odf_ol = calcDensity(ebsd('f').orientations,'halfwidth',10*degree);
odf_opx = calcDensity(ebsd('e').orientations,'halfwidth',10*degree);
odf_cpx = calcDensity(ebsd('d').orientations,'halfwidth',10*degree);

%% 7.  Define pole figures for the phases

PFs_olivine_uvw = [ ...
    Miller(1,0,0,ebsd('f').CS,'uvw'), ...
    Miller(0,1,0,ebsd('f').CS,'uvw'), ...
    Miller(0,0,1,ebsd('f').CS,'uvw'), ...
    ]

PFs_enstatite_uvw = [ ...
    Miller(1,0,0,ebsd('e').CS,'uvw'), ...
    Miller(0,1,0,ebsd('e').CS,'uvw'), ...
    Miller(0,0,1,ebsd('e').CS,'uvw'), ...
    ]

PFs_diopside_hkl = [ ...
    Miller(1,0,0,ebsd('d').CS,'hkl'), ...
    Miller(0,1,0,ebsd('d').CS,'hkl'), ...
    Miller(0,0,1,ebsd('d').CS,'hkl'), ...
    ]

%% 8. Plot pole figures olivine

figure (1)
plotPDF(odf_ol,PFs_olivine_uvw,'upper','resolution',5*degree,'contourf','colorrange','equal','minmax')
mtexColorbar

%% 9. Plot pole figures enstatite

figure (2)
plotPDF(odf_opx,PFs_enstatite_uvw,'upper','resolution',5*degree,'contourf','colorrange','equal','minmax')
mtexColorbar

%% 9.1 Plot pole figures diopside


%% 10. Compute the average stiffness tensor for each phase individually

[CVoigt_ol, CReuss_ol, CHill_ol]    = mean(C_ol,odf_ol);
[CVoigt_opx, CReuss_opx, CHill_opx] = mean(C_opx,odf_opx);
[CVoigt_cpx, CReuss_cpx, CHill_cpx] = mean(C_cpx,odf_cpx);

%% 11. Calculate the volume of each phase in this map

vol_ol  = length(ebsd('f')) ./ length(ebsd('indexed'));
vol_opx = length(ebsd('e')) ./ length(ebsd('indexed'));
vol_cpx = length(ebsd('d')) ./ length(ebsd('indexed'));

%% 12. Calculate the elastic constants for the aggregates

[CVoigt, CReuss, CHill] = mean([CVoigt_ol, CVoigt_opx, CVoigt_cpx],...
  'weights',[vol_ol, vol_opx, vol_cpx]);

%% 13. Plot the seismic velocities

plotSeismicVelocities(CHill)

%% 





 

 
 
 
 



  
  




