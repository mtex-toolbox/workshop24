%%                      Leture 1 - Basics
%
%%                      MTEX Workshop 2021
%
% This is the first script of this workshop. It should explain the basic
% principles of Matlab and MTEX.
%
%% This is the heading of the first section 
%
% We may explain first what we want to do. E.g. some Matlab basics

%% generate lists

x = 1:10;


%% pointwise operations

% generate a vector of all cube numbers
y = sin(x)


%% indexing

y([2,3,7])


%% logical indexing

isGreaterZero = y>0;

y(isGreaterZero)

%%

y(y>0)


%% change lists

y(y<0) = 0

%%
% and remove all even square numbers

y(y==0) = []


%% Lets get MTEX into play
%

% startup MTEX
run ~/mtex/master/startup_mtex.m


%%


v = vector3d.rand(10)

%%

plot(vector3d.rand(1000),'upper')

hold on
plot(vector3d.rand(1000))

plot(vector3d(1,1,1),'upper','Marker','s','MarkerSize',10,...
  'MarkerFaceColor','orange','MarkerEdgeColor','black',...
  'label','(111)','FontSize',20,'BackgroundColor','w')


hold off

%%

nextAxis
plot(vector3d.rand(1000))

%%

v = vector3d(1,1,-1,'antipodal');

plot(v,'complete')







