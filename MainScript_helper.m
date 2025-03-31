% This script was created to clean MainScript.m
% This can be run in the server to create the data. The data was uploaded
% to OSF and now it can be downloaded locally to work on analyses and plots. 
cr         = struct();
cr.codeDir = sdRP;

% WHERE THE NEW DATA IS
cr.dirs.BASE     = '/acorn/data/neuro/gari/PRF-StimDependence';
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.local    = fullfile(cr.dirs.BASE,'local');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','defineProjectDefaults');
cr.dirs.TMP      = fullfile(cr.codeDir,'DATA','ANALYSIS','TMP');
cr.dirs.FIG      = fullfile(cr.codeDir,'DATA','figures');
cr.dirs.FIGPNG  = fullfile(cr.dirs.FIG,'png');
cr.dirs.FIGSVG  = fullfile(cr.dirs.FIG,'svg');
if ~isfolder(cr.dirs.TMP); mkdir(cr.dirs.TMP); end
if ~isfolder(cr.dirs.FIG); mkdir(cr.dirs.FIG); end
if ~isfolder(cr.dirs.FIGPNG); mkdir(cr.dirs.FIGPNG); end
if ~isfolder(cr.dirs.FIGSVG); mkdir(cr.dirs.FIGSVG); end

% CONTINUE WITH THE NORMAL PROCESSING
% add to path the required matlab files inside the project, with info to run the project
% addpath(genpath(fullfile(cr.dirs.ANALYSIS,'matlabfiles')));

% Rosemary relied on this file that contains most of the subjects and other
% lists. Make it work with relative paths and store it in each project repository
% This file was used as well for:
% - copying files to a new location
% - editing mrSession to reflect the file changes
cr.bk = bookKeeping(cr);

sd_data = struct;
sd_data.CNI = struct;
sd_data.ISRAEL = struct;

sd_data.CNI.list_subInds  = 1:20;
sd_data.ISRAEL.list_subInds  = [31:36 38:44];


%% (1) Run PRFs again
% If run again, it will add a new -fFit at the end of results
% Select a single index to run just one of few subjects
% Run_PRFs(cr, sd_data.CNI.list_subInds)
% Run_PRFs(cr, sd_data.ISRAEL.list_subInds)

%% (2) READ/PREPARE DATA: WORDS (ENG/HEB), CHECKERS AND FALSEFONTS
% Generate the rmroicell that we will use in all plots in this script
% It will read results obtained by Rosemary or the re-runs in (2021,2022)
sd_data.CNI.list_roiNames = {'WangAtlas_V1d_left'
                             'WangAtlas_V2d_left'
                             'WangAtlas_V3d_left'
                             'WangAtlas_V1v_left'
                             'WangAtlas_V2v_left'
                             'WangAtlas_V3v_left'
                             'WangAtlas_hV4_left'
                             'WangAtlas_VO1_left'
                             'WangAtlas_V3A_left'
                             'WangAtlas_IPS0_left'
                             'WangAtlas_IPS1_left'
                             'WangAtlas_V1d_right'
                             'WangAtlas_V2d_right'
                             'WangAtlas_V3d_right'
                             'WangAtlas_V1v_right'
                             'WangAtlas_V2v_right'
                             'WangAtlas_V3v_right'
                             'WangAtlas_hV4_right'
                             'WangAtlas_VO1_right'
                             'WangAtlas_V3A_right'
                             'WangAtlas_IPS0_right'
                             'WangAtlas_IPS1_right'};
sd_data.ISRAEL.list_roiNames = sd_data.CNI.list_roiNames;

% Name of data type
sd_data.CNI.list_dtNames     = {'Checkers','Words','FalseFont'};
sd_data.ISRAEL.list_dtNames  = {'Checkers','Words_English','Words_Hebrew'};

% Name of results file: we have the original from Rosemary without -fFit,
% she removed it by hand. Then fFit and -fFit-fFit
sd_data.CNI.list_rmNames  = {'retModel-Checkers-css-fFit.mat'
                             'retModel-Words-css-fFit.mat'
                             'retModel-FalseFont-css-fFit.mat' };
sd_data.ISRAEL.list_rmNames  = {'retModel-Checkers-css-fFit.mat'
                                'retModel-Words_English-css-fFit.mat'
                                'retModel-Words_Hebrew-css-fFit.mat' };

% Read the data types and ROIs seen above and create a mat file. 
% Files with fits
sd_data.CNI.rmroiFname = ['rmroicell_CNI_subInds-1-20_dtNames-CB-WE-FF' ...
                          '_ROI-all_fits-fFit_2025.mat'];
sd_data.ISRAEL.rmroiFname = ['rmroicell_ISRAEL_subInds-31-36-38-44' ...
                             '_dtNames-CB-WE-WH_ROI-all_fits-fFit_2025.mat'];
% Files with time series
sd_data.CNI.TSFname = 'TS_CNI_subInds-1-20_dtNames-CB-WE-FF_ROI-all.mat';
sd_data.ISRAEL.TSFname = 'TS_ISRAEL_subInds-31-36-38-44_dtNames-CB-WE-WH_ROI-all.mat';

% Get the data (locally or online, or generate it locally)
% sd_data = read_download_data(sd_data, cr, 'rmroiCell');
% If working with time series, load them too
% sd_data = read_download_data(sd_data, cr, 'TS');
all_sd_data_file = fullfile(cr.dirs.BASE, 'PRFstim_sd_data_31march2025.mat');
% save(all_sd_data_file,'sd_data')
if ~isfile(all_sd_data_file); check_download_from_osf(fpath); end
load(all_sd_data_file)


% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames = new_list_roiNames;
% data types we want to look at
cr.defaults.covfig.vfc.list_dtNames = list_dtNames;
% names of the rm in each dt
cr.defaults.covfig.vfc.list_rmNames = list_rmNames;
% subinds = [31:36 38:44]; % Hebrew
% cr.defaults.covfig.vfc = ff_vfcDefault_Hebrew();
