hereNow = pwd;

% Add all local nested paths:
addpath(genpath(hereNow));

% Add OCDA toolbox:
theOCDApath = '~/DropboxSydneyUni/CodeToolboxes/OverlappingCommunityDetection';
addpath(theOCDApath)
cd(theOCDApath)
startup;
cd(hereNow)
