%% STEPS
% This repository is based on the code created by Rosemary Le, part of coverageReading
% repository. For the Stimulus Dependence paper (2019) we tried to make it into
% a reproducible and reusable process. The thing is that we need to separate
% this specific project from others, we need to separate data from code. 
close all; clear all;
mrvCleanWorkspace;
cr         = struct();
cr.codeDir = crRootPath;

% WHERE THE NEW DATA IS
cr.dirs.BASE     = '/black/localhome/glerma/TESTDATA/PRF-StimDependence';

% cr.dirs.BASE    = fullfile(crRootPath,'DATA');
rmpath(genpath(fullfile(crRootPath,'DATA')))
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.dirs.ANALYSIS,'matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.dirs.ANALYSIS,'matlabfiles','defineProjectDefaults');
cr.dirs.FIG      = fullfile(cr.dirs.ANALYSIS,'figures');
cr.dirs.FIGPNG  = fullfile(cr.dirs.FIG,'png');
cr.dirs.FIGSVG  = fullfile(cr.dirs.FIG,'svg');
if ~isfolder(cr.dirs.FIG); mkdir(cr.dirs.FIG); end
if ~isfolder(cr.dirs.FIGPNG); mkdir(cr.dirs.FIGPNG); end
if ~isfolder(cr.dirs.FIGSVG); mkdir(cr.dirs.FIGSVG); end

% CONTINUE WITH THE NORMAL PROCESSING
% add to path the required matlab files inside the project, with info to run the project
addpath(genpath(fullfile(cr.dirs.ANALYSIS,'matlabfiles')));

% Rosemary relied on this file that contains most of the subjects and other
% lists. Make it work with relative paths and store it in each project repository
% This file was used as well for: 
% - copying files to a new location
% - editing mrSession to reflect the file changes
cr.bk = bookKeeping(cr);

%% Run PRFs again

% subjects we want to do this for
list_subInds        = [31:36 38:44];  % Hebrew
list_subInds        = [1:20];  % Original 20
% mw (13) for Words failed, continue with the next ones for now
% list_subInds        = [18:20];
%17 and 13 failed at beginning


% Fix it: 
% list_subInds        = [13,17];

for subind =list_subInds
    
    % subind = 13;
    
    mrvCleanWorkspace;
    % subind  = list_subInds(ns);
    subname = cr.bk.list_sub{subind};
    [~,anatName]=fileparts(cr.bk.list_anatomy{subind});
    fprintf('\nSubDetails:\nInd:%i, StrInd:%s, subname:%s, Name:%s, anatName:%s\n',...
        subind,cr.bk.list_subNumberString{subind},subname,...
        cr.bk.list_names{subind},anatName)
    % Change dir, we need to run analysis where mrSession is
    chdir(cr.bk.list_sessionRet{subind})
    
    %% PRF analysis
    % Read the generic params for all subjects
    run(fullfile(cr.dirs.DEF,'prfrun_defaults.m'));
    cr.defaults.prfrun.params = params; 
    cr.defaults.prfrun.p      = p; 
    clear('params'); clear('p');
    % Read prfRun_params specific to this subject
    % run(cr.bk.list_prfParams{subind}); NOT NECESSARY
    prf.dirVistacc = cr.bk.list_sessionRet{subind};
    prf.dirAnatomy = cr.bk.list_anatomy{subind};
    prf.list_rmName= cr.bk.list_rmName{subind};
    prf.p.stimSize = cr.bk.list_stimSize(subind);
    prf.wSearch    = cr.bk.list_wSearch(subind);
    prf.prfModel   = cr.bk.list_prfModels{subind};
    prf.rois       = cr.bk.list_ROIs{subind};
    
    cr.subj.(subname).params.prf = prf;
    clear('prf');
    % This was on generics but requires specifics so... this is why I am
    % calling generics as many times as calling different subjects just in case
    cr.defaults.prfrun.params.stimSize = cr.subj.(subname).params.prf.p.stimSize; 
    % Run the prfModel with mrVista
    % RUN USING mrVISTA NORMAL INSTALLATION
        cr = cr_prfRun(cr, subind);
        % Clean workspace of globals after each subject finishes
        mrvCleanWorkspace;
    % RUN USING container prfanalyze-vista:2.0.0 (no modelpred, we get r2)
        % Generate the config file     
        % Run the container
        % pmLaunchDockerCommand('prfanalyze','ellipse','tr1dur300v3','afni6')
        % Convert the data back so that the rest of the scripts continue working
end

%% PREPARE DATA
% Generate the rmroicell that we will use in all plots in this script
% This will read the results obtained by Rosemary or the re-run in 2021
list_subInds  = [1:20];
list_roiNames = {'WangAtlas_V1v_left'
                 'WangAtlas_V2v_left'
                 'WangAtlas_V3v_left'
                 'WangAtlas_hV4_left'
                 'WangAtlas_VO1_left'
                 'lVOTRC' 
                 'WangAtlas_IPS0'
                 'WangAtlas_IPS1'};
list_dtNames  = {'Checkers','Words','FalseFont'};
list_rmNames  = {'retModel-Checkers-css-fFit.mat'
                 'retModel-Words-css-fFit.mat' 
                 'retModel-FalseFont-css-fFit.mat' };
% {
% Use the originals calculated by Rosemary
list_rmNames  = {'retModel-Checkers-css.mat'
                 'retModel-Words-css.mat' 
                 'retModel-FalseFont-css.mat'};
%}
list_rmDescripts = {'Words'...  % Words (large bars)
                    'Checkers'...
                    ... % 'Words_English'...
                    ... % 'Words_Hebrew'... % Words (smalls bars)
                    'FalseFont'};

rmroiCell=ff_rmroiCell(cr,list_subInds,list_roiNames,list_dtNames,list_rmNames);
% Save rmroicell just in case
% save(fullfile(crRootPath,'DATA',...
%      'rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-Rosemary.mat'),'rmroiCell')
% load(fullfile(crRootPath,'DATA',...
%      'rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-Rosemary.mat'),'rmroiCell')

%% FIGURE 1: (C) Groups coverage plots

% (A) Explain how to obtain


% (B) Explain how to obtain


% 
% With the new data the groups plots look different, but it seems that it
% is due to thresholds
% >> Check colormap limits  so that checker looks bigger than words

% Group COVERAGE plots, take all subjects from list_subInds

% Select subjects we want to plot
% subinds = [1:20]; % Stanford Subjects, 1 is gomez, find anatomicals
% subinds = [1:12,14:16,18:20];

% Read the generic params for coverage for all subjects
cr.defaults.covfig.vfc = ff_vfcDefault();
cr.defaults.covfig.vfc.list_roiNames = list_roiNames;
% data types we want to look at
cr.defaults.covfig.vfc.list_dtNames = list_dtNames;
% names of the rm in each dt
cr.defaults.covfig.vfc.list_rmNames = list_rmNames;
% subinds = [31:36 38:44]; % Hebrew
% cr.defaults.covfig.vfc = ff_vfcDefault_Hebrew();

% Launch the function
figFunction_coverage_maxProfile_group(cr,list_subInds,'flip',false, ...
                                      'savefig',true, 'vers','v01_oldfit',...
                                      'rmroiCell',rmroiCell)
                               
%% FIGURE 2: (C) Scatterplots: word-checkerboard
rmroiCell_WC    = rmroiCell(:,1:6,1:2);
list_roiNames16 = list_roiNames(1:6);
fname           = 'scatterplot_WordVsCheck_6ROIs_20subs_RosemaryFit_v01';
[rmroiCellSameVox,C_data,cmapValuesHist,maxValue]=crThreshGetSameVoxel(...
                                                                    cr,...
                                                          rmroiCell_WC,...
                                                          list_subInds,...
                                                       list_roiNames16,...
                                                          list_rmNames,...
                                                      list_rmDescripts,...
                                                                   fname);

%% FIGURE 3: (B) Line plots


%% COVERAGE: individual plots 
for subind = 1 [1:12,14:16,18:20] % list_subInds
    subname = cr.bk.list_sub{subind}
    %% Plot the coverage figures
    % Read the coverage figure params
    % run(cr.bk.list_coverageFigure_defaults{subind});
    % Defaults
    covfig.vfc              = ff_vfcDefault;
    covfig.titleDescript    = 'FOV';
    % vfc threshold
    covfig.cothresh         = 0.2; 
    covfig.vfc.cmap         = 'hot';
    covfig.vfc.addCenters   = true;
    covfig.vfc.contourPlot  = true;
    % ROIs
    covfig.list_roiNames    = {'lVOTRC'};
    covfig.list_roiNames    = {'WangAtlas_VO1_left'};
    % dt and rm names
    covfig.list_dtNames     = {'Checkers','Words','FalseFont'};
    covfig.list_rmNames     = {'retModel-Checkers-css-fFit.mat'
                               'retModel-Words-css-fFit.mat'
                               'retModel-FalseFont-css-fFit.mat'};
    % {
    % Old original fits, basically they are the same
    covfig.list_rmNames     = {'retModel-Checkers-css.mat'
                               'retModel-Words-css.mat'
                               'retModel-FalseFont-css.mat'}; 
    %}
    covfig.list_rmDescripts = {'Checkers', 'Words','FalseFont'};
    
    
    
    cr.subj.(subname).params.covfig = covfig;
    clear('covfig');
    % Plot it
    [RFcov,weight, data] = figFunction_coverage_individual(cr, subind);
    % Clean workspace of globals after each subject finishes
    mrvCleanWorkspace;
end



%% Notes
% + this is lVOTRC, do the same for V1-4,hvo1 (the ones in the paper)
% - two plots, or separate long versus short lines
% + go to white background
% - obtain numbers that show that eccc>eccw is basically noise (it will be
%      only a caption in the figure)
% + DO NOT plot any ecc diff of +- 0.5 deg
% - In the scatterplot with the light blue cones:
%    + Remove cone
%    + add +-1.5deg band
%    + below, we will only have about 3% of the data
% - compare the - and + differences, they shuold be the same for V1 v2 and
% then start changing
% - when re-running the fits, use two different HRFs (simulate that they
% will actually have size differences) and show that the effect and the
% centers will not vary
% - test: select and HRF that gives the correct size in V3 at 5deg eccen,
% and run all the analyses with this HRF
% - Use all V1, not only ventral


% test of radiality from 5 to 7 deg and 8 to 12 degs, move CB to the horizontaal and maintain the angle for words



% TODO: 
% - print the line plots again  and finish figure 3
% - plot the results in IPS to check if they hold

%% DATA ANALYSIS DONE BY ROSEMARY
%{

%% (once) Upload data
% Original data in black.stanford.edu
% Create local structure that we couls upload to FW later
% - Prepare using the /utilities/prepareDataUploadFW.py script
% - Upload to FW using the fw command line utility
%      fw import folder .
%   Leonardo did it from his computer
% - TODO: add the subject and acqu metadata

%% (once) Run the analysis and upload the data table to the collection
% This example below is what we use for DWI analysis with RTP. Now we need
% something similar for PRF-s: run them in the server, and download only the
% relevant data for our analyses and plots (and next we could obtain the
% analyses online as well. 

%   Usually I run this in Google server and then continue locally
%{
% Run the analyses
serverName     = 'stanfordlabs';
collectionName = 'HCP_Depression';
gearName       = 'afq-pipeline'; 
gearVersion    = '3.0.7';
% Before launching check the script for the correct analysis config params
dr_fwLaunchJobs(serverName, collectionName, gearName, gearVersion)

% Create the datatable (only fa for now)
measurement    = 'fa';
% get all the analysis in the collection
JL = dr_fwCheckJobs(serverName, collectionName);
% filter the analyses
state         = 'complete';
dateFrom      = '04-Feb-2019 00:00:00';
labelContains = 'v3.0.7:';
t  = JL(JL.state==state & JL.gearName==gearName & ...
       JL.gearVersion==gearVersion & JL.JobCreated>dateFrom & ...
       contains(string(JL.label), labelContains),:);
% read data of interest and create table
dt = dr_fwReadDtFromAnalysisTable(serverName, t, measurement);
% It takes a lot of time, save it locally...
localfname    = fullfile(stRootPath,'local','tmp', ...
                      sprintf('AllV02_HCP_Depression_%s.mat',measurement));
save(localfname, 'dt')
%  ...and upload it to the collection
st   = scitran(serverName); st.verify;
cc   = st.search('collection','collection label exact',collectionName);
stts = st.fileUpload(localfname, cc{1}.collection.id, 'collection');
% Check that the data is there
[~,fname,ext] = fileparts(localfname);
data          = load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
                                                  [fname ext],localfname));
%}

%% Download the data from the collection for analysis
%{
% Every time we want to re-run the data, we can download it from the server. 

% Download or clone the repository
%       !git clone https://github.com/garikoitz/paper-HCPDEPRESSION.git
% 
clear all; close all; clc;
% Add the root of the repository to the Matlab path (or run this code):

    % cd('<path-to-your-code>/paper-reproducibility')
    cd('~/soft/paper-HCPDEPRESSION')
    rootDir = pwd;
    addpath(genpath(rootDir));

% Specify a path to save the output figures. 
paperPath = '~/gDrive/STANFORD/PROJECTS/2019 Depression and WM (Leonardo-Gari-Brian)';
saveItHere = string(fullfile(paperPath, 'VERSION_01/figures/sources'));


% Read the data 
% Check if there is a local cache, otherwise download it from FW
DataVersion    = '01';
collectionName = 'HCP_Depression';
measure        = 'fa';

fname          = sprintf('AllV%s_%s_%s.mat',DataVersion, collectionName, measure);
localfname = fullfile(paperReprPath,'local',fname);
if exist(localfname,'file')
    data = load(localfname);
else  % Download it from the Flywheel collection attachment
    serverName     = 'stanfordlabs';
    st  = scitran(serverName);
    cc  = st.search('collection','collection label contains',collectionName);
    data=load(st.fw.downloadFileFromCollection(cc{1}.collection.id,fname,localfname));
end

% This script can be run directly using the Run button or step by step. 
%}

%% Data preparation
if (0)
% Do not do it, use what we have already    
    
% apply canonical x form -- for every nifti
niftiApplyCannonicalXform

% acpc align the anatomical 
mrAnatAverageAcpcNifti

% run freesurfer
eval(['! recon-all -i ' pathT1 ' -subjid ' dirNameFreesurfer ' -all'])

% ribbon from freesurfer into class file -- t1_class.nii.gz
fs_ribbon2itk(inputRibbonFile, outputClassNii, [], pathT1, [])

% TODO: see if we can substitute the previous steps using fmriprep
end

%% Initialization 
if (0)
% Do not do it, use what we have already    
   
% Initialize. This will create the mrSESSION in the root folder. 
cr_mrInit(cr, opt);
end

%% Allignment and dataType creation
if (0)
% Do not do it, use what we have already    
   
% align the inplane to anatomical
s_alignInplaneToAnatomical;

% This is the message when I closed it:
%    The alignvolumedata GUI has been closed. in case you forgot to export the transformation, here it is:
%    tr = maketransformation([0 0 0],[1 2 3],[95.9343856857804 106.096638132966 92.1256300034666],[1 2 3],[-88.7924000052859 -2.12048422951174 98.6511462809804],[192 192 62],[208.000007629395 208.000007629395 62],[0.998614818358445 -1.00108310465201 1.9968769575349],[0 0 0],[0.000745334952393827 -0.000334517493656755 0.000639298200901277],[0 0 0]);
% make the transformation into a 4x4 matrix
%    rx.volVoxelSize = [2 2 2];
%    T = transformationtomatrix(tr,0,rx.volVoxelSize);
% vw = initHiddenInplane; mrGlobals; 
% mrSESSION.alignment = T;
% saveSession;




% TODO??
% specify segmentation file, go to gray view to run the prfs

% create a new dataTYPE which is the average of the 4 runs (note, still in INPLANE)
% this script will also xform the data into the gray

% the names of the dataTYPES we want to create
dtsToCreate = {
    'WordFalse1';     % 1 % the 1st checker and the 2nd word
    'WordFalse2';     % 2 % the 2nd checker and the 1st word
    };
opt.dtsToCreate = {
    'Checkers1'         % 1
    'Words_English1'    % 2  
    'Words_Hebrew1'     % 3
    'Words_Hebrew2'     % 4
    };
% 
% The datatype the scan belongs to. For example, a 1 means that the first
% scan is in the first dataTYPE specified in dtsToCreate
% use 0 if the scan is not used in the creation of a new dt
opt.dtAssignments = [
    0;
    0;
    1;
    2;
    3;
    4;
    ];

% make the new tseries from the most processed time series
opt.dtToAverage = 'MotionComp_RefScan1';
% Run it
cr = crNewDataTypes(cr, opt); 
end

%% Stimuli preparation
if (0)
% Do not do it, use what we have already    
   
% make a Stimuli folder in the same place as the mrSESSION.mat
% for localizer GLM analyses, make Stimuli/Parfiles
% 2 things go into Stimuli
% params file -- mrVista writes this to desktop
% image matrix -- (part of the params file)
end

%% Analysis of the data
if (0)
% The code here will go to a Docker container. Make the figures in the paper reproducible
% The code below will come as a combination from pmVistasoft.m that I did
% adapting a script from Jon Winawer, and s_prfRun.m Rosemary Le's script. 

cr.dirs.prfRun_params()



% subjects we want to do this for
opt.list_subInds = [31]; 
% dataTYPE name. Can run for mutiple datatypes 
opt.list_rmName = {'Words_Hebrew1','Words_English1'}; 
% roi name. assumes in shared anatomy directory (change this to be self contained)
% if we want to run on the whole brain, assign this the empty string ''
% assign this to be a string in a cell otherwise {'LV1_rl'}
opt.list_rois = {'lVOTRC'}; 
% prf model. Specify in a cell. Options: 
% {'one oval gaussian' | 'onegaussian' | 'css'}
% Note: if we want to specify multiple models, change the naming
% convention. See outFileName
opt.prfModel = {'css'}; 
% search type. 
% 1 = grid search only ("coarse"),
% 2 = minimization search only ("fine"),
% 3 = grid followed by minimization search [default]
%   note, there is another option which is to find the hrf as well
opt.wSearch = 3; 
% radius of circle retinotopy in visual angle degrees
opt.p.stimSize = 7; %%% Is this true for the Hebrew subject we just selected?
% define things common to all datatypes
% name of params file
opt.p.paramsFile_Knk        = 'Stimuli/params_knkfull_multibar_blank.mat';  % Words and FalseFont
opt.p.paramsFile_Checkers   = 'Stimuli/params_checkers.mat';                % Checkers
% image file
opt.p.imFile_Knk            = 'Stimuli/images_knk_fliplr.mat';              % Words and FalseFont
opt.p.imFile_Checkers       = 'Stimuli/images_8barswithblank_fliplr.mat';   % Checkers
% params common to all dts
opt.params.stimSize         = opt.p.stimSize; 
opt.params.fliprotate       = [0 0 0]; 
opt.params.stimType         = 'StimFromScan';
opt.params.stimWidth        = 90;               
opt.params.stimStart        = 0;                
opt.params.stimDir          = 0;                
opt.params.nCycles          = 1;               
opt.params.nStimOnOff       = 0;                
opt.params.nUniqueRep       = 1;                
opt.params.nDCT             = 1;     
opt.params.hrfType          = 'two gammas (SPM style)';
opt.params.hrfParams        = {[1.6800 3 2.0500] [5.4000 5.2000 10.8000 7.3500 0.3500]}; 
opt.params.imfilter         = 'binary';
opt.params.jitterFile       = 'Stimuli/none';

% Run it

% We downloaded the mrSESSION from black, edit it before using it locally
results = cr_prfRun(cr, opt);
end

%% Generate the individual figure
if (0)
% {
% For a full list of options, see edit IndivFigDefaults
% title description
opt.titleDescript    = 'FOV';
% vfc threshold
opt.vfc              = ff_vfcDefault_Hebrew; 
opt.vfc.cmap         = 'hot';
opt.vfc.addCenters   = true;
opt.vfc.contourPlot  = true;
% subjects
opt.list_subInds     =  [31];
% session
% opt.list_path        = cr.bk.list_sessionRet; 
% opt.list_path        = list_sessionRet; 
% ROIs
opt.list_roiNames    = {'lVOTRC'};
% dt and rm names
opt.list_dtNames     = {'Words_Hebrew','Words_English'};
opt.list_rmNames     = {'retModel-Words_Hebrew-css.mat'
                        'retModel-Words_English-css.mat'};
opt.list_rmDescripts = {'Words_Hebrew', 'Words_English'};
figFunction_coverage_individual(cr, opt);
%}
end

%}

%% Line plots: from checquerboard to word
% Combine code in figScript_centers_* to create the plots Brian nd I
% discussed
%   - left VOTC, line plot going from word to cb
%   - black if both word and cb are in the same quadrant (); else gray
%   - continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
%   - count percentage of continuous black lines over rest:
%      - bin by size (small sizes very unreliable)
%      - 

% SELECT SUBJECTS AND MODELS
%{
onlyStanford  = [1:20];  % Why not the rest? Ask MBS/RL
onlyHebrew    = [31:36 38:44]; 
list_subInds  = [onlyStanford onlyHebrew]; 
list_subInds  = [1:12]; 
list_subInds = [1:12,14:16,18:20];

list_path     = cr.bk.list_sessionRet; 
list_roiNames = {'WangAtlas_V1v_left';
                 'WangAtlas_V2v_left';
                 'WangAtlas_V3v_left';
                 'WangAtlas_hV4_left';
                 'WangAtlas_VO1_left';
                 'lVOTRC';
                 'WangAtlas_IPS0';
                 'WangAtlas_IPS1'};
% list_roiNames = {'LV1_rl'
%                  'LV2v_rl'
%                  'LV3v_rl'
%                  'LhV4_rl'
%                  'LVO1_rl'
%                  'lVOTRC' };
list_dtNames = {'Words'...
                'Checkers'...
                'Words_English'...
                'Words_Hebrew'...
                'FalseFont'};
% ret model names
list_rmNames = {'retModel-Words-css-fFit.mat'...
                'retModel-Checkers-css-fFit.mat'...
                'retModel-Words_English-css.mat'...
                'retModel-Words_Hebrew-css.mat'...
                'retModel-FalseFont-css.mat'};
list_rmDescripts = {'Words'...  % Words (large bars)
                    'Checkers'...
                    'Words_English'...
                    'Words_Hebrew'... % Words (smalls bars)
                    'FalseFont'};
%}









