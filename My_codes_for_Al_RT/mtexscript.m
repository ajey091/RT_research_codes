%% Import Script for EBSD Data
%
% This script was automatically created by the import wizard. You should
% run the whoole script or parts of it in order to import your data. There
% is no problem in making any changes to this script.

%% Specify Crystal and Specimen Symmetries
clear
close all
tic
% crystal symmetry
CS = crystalSymmetry('Fm-3m','mineral','Au');   

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Specify File Names

pname = '/Users/b119user/Downloads/Ajey/Al_Ajey';


% which files to be imported
fname = '1.ang';

%% Import the Data

% create an EBSD variable containing the data
% ebsd = loadEBSD(fname,'CS',CS,...
%                 'interface','generic','Bunge','ignorePhase',[0 2],...
%                  'ColumnNames', { 'Phase' 'x' 'y' 'Euler 1' 'Euler 2' 'Euler 3'},...
%                  'Columns', [2 3 4 5 6 7])
ebsd = loadEBSD(fname,CS,'interface','ang',...
  'convertSpatial2EulerReferenceFrame');
% close all


% saveas(gcf,'Au_IPF.png');
grains = calcGrains(ebsd,'angle',5*degree);
oM = ipdfHSVOrientationMapping(ebsd);
%oM.inversePoleFigureDirection = mean(ebsd.orientations) * oM.whiteCenter;
oM.inversePoleFigureDirection = zvector;
oM.colorStretching = 1;
setMTEXpref('xAxisDirection','north')
setMTEXpref('zAxisDirection','intoPlane')
setMTEXpref('showMicronBar','off')
% plot(oM,'FontSize',32)
% 
% hold on
% % plotIPDF(ebsd('Au').orientations,xvector,'markerSize',2,'points',20000,'marker','o','markerfacecolor','none','markeredgecolor','k')
% plotIPDF(grains('Au').meanOrientation,xvector,'markerSize',3,'points',400,'marker','o','markerfacecolor','none','markeredgecolor','k')

% define the meanFilter
F = medianFilter;

% define the size of the window to be used for finding the median
F.numNeighbours = 3; % this corresponds to a 7x7 window
% smooth the data
ebsd_smoothed = smooth(ebsd,F);

grains_smoothed = calcGrains(ebsd_smoothed,'angle',5*degree);
% plot(oM,'FontSize',32)
% 
% hold on
% % plotIPDF(ebsd('Au').orientations,xvector,'markerSize',2,'points',20000,'marker','o','markerfacecolor','none','markeredgecolor','k')
% plotIPDF(grains_smoothed('Au').meanOrientation,xvector,'markerSize',3,'marker','o','markerfacecolor','none','markeredgecolor','k')

figure
plot(ebsd,oM.orientation2color(ebsd.orientations))

figure
plot(ebsd_smoothed('Au'),oM.orientation2color(ebsd_smoothed('Au').orientations))
% saveas(gcf,'EBSD.png')
% Getting some control over what's going on
% smoothed_x = ebsd_smoothed('Au').x;
% smoothed_y = ebsd_smoothed('Au').y;
% smoothed_ori = ebsd_smoothed('Au').orientations.Euler;
% IPFcolors = oM.orientation2color(ebsd_smoothed('Au').orientations);
% figure('pos',[10,10,1200,600])
% scatter(smoothed_y,smoothed_x,[],IPFcolors)


toc