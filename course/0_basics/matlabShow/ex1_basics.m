%%                        MTEX Workshop 2021
%
% This is the first script of this workshop. It should explain the basic
% principles of Matlab and MTEX. 
%
%% This is the heading of the first section 
%
% We may explain first what we want to do. E.g. some Matlab basics


%% List generation

% generate a vector of all square numbers up to 25
x = (1:5).^2

%% column vectors vs. row vectors

x.'

%% compute with vectors

x .* x

%%

x * x.'


%% Sublists

% select by index
x(3)

% select by condition
isEven = round(x/2) == x/2

x(isEven)

%% List manipulation

% change entries


% combine lists

% remove entries
x(iseven(x)) = []


%% Lets come MTEX into play
%

% startup MTEX
run ~/mtex/master/startup_mtex.m


%% vector3d


%% lists of vectors


v = vector3d.X


%% import vector data from file

% specify filename
fname = 'smiley.txt';

% specify column names in file
cNames = {'polar','azimuth','value'};

% import data from file
[v,data] = vector3d.load(fname,'columnNames',cNames);

%% plot data

plotx2north
plot(v,data.value,'upper')

% rotate data
%plot(rotate(v,30*degree),data.value,'upper')

%%

mtexColorMap blue2red
mtexColorbar

%% export data to file
export(v,'test.txt',data)

