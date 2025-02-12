%% 
% Lets consider a quartz crystal. 

cs = crystalSymmetry.load('quartz')

%%
% Its shape is mainly bounded by the following faces

m = Miller({1,0,-1,0},cs);  % hexagonal prism
r = Miller({1,0,-1,1},cs);  % positive rhomboedron, usally bigger then z
z = Miller({0,1,-1,1},cs);  % negative rhomboedron
s1 = Miller({2,-1,-1,1},cs);% left tridiagonal bipyramid
s2 = Miller({1,1,-2,1},cs); % right tridiagonal bipyramid
x1 = Miller({6,-1,-5,1},cs);% left positive Trapezohedron
x2 = Miller({5,1,-6,1},cs); % right positive Trapezohedron

%%

n = [m,r,z,s1,s2,x1,x2];

for k = 1:length(n)

  plot(n(k).symmetrise,'projection','eangle','upper',...
    'displayName',char(n(k)))
  %plot(n(k).symmetrise,'3d', 'displayName',char(n(k)))
  hold on

end

hold off

legend('show','Location','eastoutside','fontSize',20)


%%
% If we take only the first three faces prism + rhomboedron
% we end up with

N = [2*m,r,z];
cS = crystalShape(N)

plot(cS,'colored')

%%

N = [2*m,r,z];

cS = crystalShape(N);
plot(cS,'colored')

%%

% collect the face normal with the right scalling
% make the positive rhomboedron bigger then the negative 
N = [2*m,r,0.9*z];

cS = crystalShape(N);
plot(cS,'colored')

%%
% Finaly, we add the tridiagonal bipyramid and the positive Trapezohedron

% collect the face normal with the right scalling
N = [2*m,r,0.9*z,0.7*s1,0.3*x1];

cS = crystalShape(N);
plot(cS,'colored')

%% Marking crystal faces
% We may colorize the faces according to their lattice planes using the
% command

cS = crystalShape.topaz
plot(cS,'colored')



%%
% or even label the faces directly

plot(cS)
N = unique(cS.N.symmetrise,'noSymmetry','stable');
fC = cS.faceCenter;

for i = 1:length(N)
  text3(fC(i),char(round(N(i)),'latex'),...
    'scaling',1.1,'interpreter','latex')
end


%% Defining complicated crystals more simple
% We see that defining a complicated crystal shape is a tedious work. To
% this end MTEX allows to model the shape with a habitus and a extension
% parameter. This approach has been developed by J. Enderlein in
% <https://library.wolfram.com/infocenter/Articles/3279 A package for
% displaying crystal morphology. Mathematica Journal, 7(1), 1997>. The two
% parameters are used to model the distance of a phase from the origin.
% Setting all parameters to one we obtain

% take the face normals unscaled
N = [m,r,z,s2,x2];

habitus = 1;
extension = [1 1 1];
cS = crystalShape(N,habitus,extension);
plot(cS,'colored')


%%
% The scale parameter models the inverse extension of the crystal in each
% dimension. In order to make the crystal a bit longer and the negative
% rhomboedrons smaller we could do

extension = [0.9 1.1 1];
cS = crystalShape(N,habitus,extension);
plot(cS,'colored')

%%
% Next the habitus parameter describes how close faces with mixed hkl are
% to the origin. If we increase the habitus parameter the trapezohedron and
% the bipyramid become more and more dominant

habitus = 1.1;
cS = crystalShape(N,habitus,extension);
plot(cS,'colored'), snapnow

figure

habitus = 1.2;
cS = crystalShape(N,habitus,extension);
plot(cS,'colored'), snapnow

figure

habitus = 1.3;
cS = crystalShape(N,habitus,extension);
plot(cS,'colored')

%% Select faces
% A specific face of the crystal shape may be selected by its normal vector

plot(cS)
hold on
plot(cS(Miller(0,-1,1,0,cs).symmetrise),'FaceColor','DarkRed') 
hold off

%% Gallery of hardcoded crystal shapes

plot(crystalShape.olivine,'colored')

%%

plot(crystalShape.garnet,'colored')

%%

plot(crystalShape.topaz,'colored')

%%

plot(crystalShape.plagioclase,'colored')
