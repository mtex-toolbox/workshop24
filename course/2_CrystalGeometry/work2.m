rot = rotation.rand(400)

%%

plot(rot,'filled')

%saveFigure('pic/rotEuler.png','crop')

%%

plot(rot,'axisAngle','linewidth',2)
axis off

% saveFigure('pic/rotAxisAngle.png')

%%


plotSection(rot,'axisAngle',[50,75,100,125,150,175].*degree)

% saveFigure('pic/rotAxisAngleSec.png')

%%

plotSection(rot,'phi2',[50,75,100,125,150,175].*degree)

% saveFigure('pic/rotEulerSec.png')

%%

Z   = vector3d.Z
rot = rotation.byAxisAngle(Z, 120*degree)
cs  = crystalSymmetry.byElements(rot)

cs.rot

%%

m  = reflection(vector3d.X)
cs = crystalSymmetry.byElements([rot, m])

cs.rot

%%

cs.add(rotation.inversion)
