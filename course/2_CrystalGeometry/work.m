
%%
% define an orthogonal coordinate system

u1 = vector3d(1,2,1); u1 = u1.normalize;
u2 = vector3d(0,-1,2); u2 = u2.normalize;
u3 = cross(u1,u2); u3 = u3.normalize;

U = matrix(u1,u2,u3)

%%

v1 = vector3d(1,1,1); v1 = v1.normalize;
v2 = vector3d(1,-1,0); v2 = v2.normalize;
v3 = cross(v1,v2); v3 = v3.normalize;

V = matrix(v1,v2,v3)

%% active rotation u -> v

rotActive = rotation.map(u1,v1,u2,v2)

matrix(rotActive)

V * U'


matrix(rotActive) * matrix(u1,u2,u3)


%% passive rotation
% translate coefficients with respect to u into coefficients with respect
% to v

% coefficient vector with respect to u
cu = [1;2;3];
w = [u1,u2,u3] * cu 

% coefficient vector with respect to v
%cv = matrix(inv(rotActive)) * cu;

cv = V' * U * cu;

w = [v1,v2,v3] * cv 

%%

rotActive = rotation.byMatrix(V * U')

rotPassive = rotation.byMatrix(V' * U)

%
% V * rotPassive' * V' =  V * (V'*U)' * V' = V * U' = rotActive
% U * rotPassive  * U' =  U * (V'*U) * U' = U * V' = (V * U')' = inv(rotActive)


%%

rotActive.angle
rotPassive.angle

%%
%
% V'*U * cu = cu -> U * cu = V * cu
% V * U' * w = w -> U' * w = V' * w ---> w = U * U' * w = U * V' * w 


rotActive.axis
rotation.byMatrix(U) * rotPassive.axis

%%

rotActive.axis
rotation.byMatrix(V) * rotPassive.axis

%% 

rot1 = rotation.rand
rot2 = rotation.rand

matrix(rot1 * rot2)
matrix(rot1) * matrix(rot2)

rotation.byMatrix(matrix(rot1))


%%

v = vector3d(1,1,1);
omega = linspace(0,2*pi);
rot = rotation.byAxisAngle(v,omega)

%%

line(rot * vector3d(1,-1,-3:3),'upper')

%%

r1 = rotation.byEuler(10*degree,20*degree,30*degree)

r2 = rotation.byAxisAngle(zvector,10*degree) * ...
  rotation.byAxisAngle(xvector,20*degree) * ...
  rotation.byAxisAngle(zvector,30*degree)


















