%% (0) INIT
% This repository is based on the code created by Rosemary Le, part of coverageReading
% repository. For the Stimulus Dependence paper (2019) we tried to make it into
% a reproducible and reusable process. The thing is that we need to separate
% this specific project from others, we need to separate data from code.
tbUse PRF-StimDependence;

close all; clear all;
mrvCleanWorkspace;
cr         = struct();
cr.codeDir = sdRP;

% WHERE THE NEW DATA IS
cr.dirs.BASE     = '/black/localhome/glerma/TESTDATA/PRF-StimDependence';
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','defineProjectDefaults');
% cr.dirs.FIG      = fullfile(cr.codeDir,'DATA','figures');
cr.dirs.FIG     = fullfile('/Users/glerma/Library/CloudStorage/GoogleDrive-garikoitz@gmail.com/My Drive/STANFORD/PROJECTS/2018 Reading across maps (Rosemary)/__PUBLISH__/2022_PNAS(3rd)','figures');
cr.dirs.FIGPNG  = fullfile(cr.dirs.FIG,'png');
cr.dirs.FIGSVG  = fullfile(cr.dirs.FIG,'svg');

cr.bk = bookKeeping(cr);

%% Create all the rmroi_Cell that we are interested
% CNI
fnameCNI = fullfile(sdRP,'DATA','rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-new_LeftRightROIs_2023.mat');
CNI = load(fnameCNI);
size(CNI.rmroiCell)
% Order of elements in the rmroiCell is: 
list_dtNames  = {'Checkers','Words','FalseFont'};

% ISRAEL
fnameHEB = fullfile(sdRP,'DATA','rmroicell_subInds-31to36-38to44_dtNames-ALL-LeftRight_fits-new_2023.mat');
HEB = load(fnameHEB);
size(HEB.rmroiCell)
% Order of elements in the rmroiCell is: 
list_dtNames     = {'Words_English','Words_Hebrew','Checkers'};


