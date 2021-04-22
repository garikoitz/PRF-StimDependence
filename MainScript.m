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
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.ORG      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','organization');
cr.dirs.DEF      = fullfile(cr.codeDir,'DATA','ANALYSIS','matlabfiles','defineProjectDefaults');
cr.dirs.FIG      = fullfile(cr.codeDir,'DATA','figures');
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
if 0
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
readExisting = true;
if readExisting
    load(fullfile(crRootPath,'DATA',...
      'rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-Rosemary.mat'),'rmroiCell');
else
    rmroiCell=ff_rmroiCell(cr,list_subInds,list_roiNames,list_dtNames,list_rmNames);
    % Save rmroicell just in case
    save(fullfile(crRootPath,'DATA',...
          'rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-Rosemary.mat'),'rmroiCell')
end

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
fname = 'Fig1_'; % '' for not saving
figFunction_coverage_maxProfile_group(cr,list_subInds, 'flip',false, ...
                                      'rmroiCell',rmroiCell,...
                                      'fname', fname, 'vers','v01_oldfit');
                               
%% FIGURE 2: (C) Scatterplots: word-checkerboard
% Order is CB, W, FF, invert it so that it is W then CB
rmroiCell_WC    = rmroiCell(:,1:6,1:2);
rmroiCell_WC    = flipdim(rmroiCell_WC,3);
list_roiNames16 = list_roiNames(1:6);

% Obtain equally thresholded voxels to scatterplot
[R,C_data,cr]=crThreshGetSameVoxel(cr,...
                                   rmroiCell_WC,...
                                   list_subInds,...
                                   list_roiNames16,...
                                   'cothres', 0.2,...
                                   'fieldrange', 15);

% Plot it
fname = 'scatterplot_WordVsCheck_6ROIs_20subs_RosemaryFit_v01';
crCreateScatterplot(R,C_data,cr,...
                    list_subInds,...
                    list_roiNames16,...
                    list_rmNames,...
                    list_rmDescripts,...
                    'ecc', ...  % 'co'
                    fname);

                
%% FIGURE 3: (B) Line plots
% Uses the same voxel calculations from the previous plot
% If doubt or this is moved, calculate it again here

% Plot it
fname = 'lineplot_WordVsCheck_6ROIs_20subs_RosemaryFit_v01';
crCreateLinePlot(R,C_data,cr,...
                 list_subInds,...
                 list_roiNames16,...
                 list_rmNames,...
                 list_rmDescripts,...
                 fname);
             




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
% end

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

%% Get the cell of rms so that we can threshold
if 0
    rmroiCell = ff_rmroiCell(cr,list_subInds, list_roiNames, list_dtNames, ...
            list_rmNames, 'list_path', list_path);
    % SAVE THIS TO WORK LOCALLY
    mkdir(fullfile(crRootPath,'DATA'))
    % 20 subs, words and CBs
    % save(fullfile(crRootPath,'DATA','rmroicell_1to20.mat'),'rmroiCell');
    % ALL subs, words and CBs and FF and Heb
    save(fullfile(crRootPath,'local','rmroicell_allandall.mat'),'rmroiCell');
end
if 0
    % load(fullfile(crRootPath,'DATA','rmroicell_1to20.mat'));
    load(fullfile(crRootPath,'local','rmroicell_allandall.mat'));
    % Only the 20
    list_subInds  = [1:20]; 
    rmroiCell = rmroiCell(list_subInds,(1:6),[1,2]);
    % Only the ones with FF
    % list_subInds  = [1:12]; 
    % rmroiCell = rmroiCell(list_subInds,(1:6),[1,5]);
    % Only the Hebrew
    % rmroiCell = rmroiCell((21:end),(1:6),[1,5]);
end
                
%%  Threshold and get identical voxels for each subject
% values to threshold the RM struct by
% vfc = ff_vfcDefault_Hebrew;
% The defaults are different for the two projects, sofor now use the most
% restrictive one
% vfc threshold
vfc.prf_size         = true;
vfc.fieldRange       = 15; % 7;
vfc.method           = 'max';
vfc.newfig           = true;
vfc.nboot            = 50;
vfc.normalizeRange   = true;
vfc.smoothSigma      = false;
% Thresholds
vfc.sigmaEffthresh   = [.2 15];  % [0  7]; % sigma effect (sigmaMajor/sqrt(exponent))
vfc.sigmaMajthresh   = [0.5 8];   % [0 14]; % sigma major (before the exponent)
vfc.cothresh         = 0.20; % Variance explained, 20%
vfc.cothreshceil     = 1; % 0.2 looking for noise. 1 for normal thresholding. don't get voxels higher than this
vfc.threshByCoh      = true;
vfc.eccthresh        = [0 vfc.fieldRange];
vfc.quadthresh       = [1 4];  % If center not on these quads, remove
% Fig control
vfc.nSamples         = 128;
vfc.meanThresh       = 0;
vfc.weight           = 'fixed';
vfc.weightBeta       = false;
vfc.cmap             = 'jet';
vfc.clipn            = 'fixed';
vfc.addCenters       = false;
vfc.verbose          = prefsVerboseCheck;
vfc.dualVEthresh     = 0;
vfc.ellipsePlot      = false;
vfc.ellipseLevel     = 0.5;
vfc.ellipseColor     = [1 0 0];
vfc.contourPlot      = true;
vfc.contourLevel     = 0.5;
vfc.contourColor     = [0 0 0];
vfc.tickLabel        = false;
vfc.contourBootstrap = false;
vfc.gridColor        = [.6 .6 .6];
vfc.backgroundColor  = [1 1 1];
vfc.cmapRange        = [0 pi]; % the range over which color bar
vfc.cmapValues       = flipud(jetCmap(0,128)); % the colorbar values
vfc.alphaValue       = ''; % 0.5;
vfc.alphaValueDot    = ''; % 0.8;
vfc.lineWidth        = 1.5; % thicker lines for transparent. 1 works well for opaque
% by definition, when eccentricity does not shift a lot, theta will not
% shift either. We can look at theta shifts in the voxels whose
% eccentricity have shifted a lot
% thetaShiftByEccThresh = true;
% eccThresh = 3;
% Update the cr structure that I created that go into the new functions
cr.defaults.covfig.vfc = vfc;
% INITIALIZE SOME THINGS
numRois = length(list_roiNames);
numSubs = length(list_subInds);
% cell for linearizing the data (a vector for each ROI)
L_data = cell(1, numRois);
X_rm1 = cell(1, numRois);
Y_rm1 = cell(1, numRois);
X_rm2 = cell(1, numRois);
Y_rm2 = cell(1, numRois);
Ecc_rm1 = cell(1, numRois);
Ecc_rm2 = cell(1, numRois);
rmDescript1 = list_rmDescripts{1};
rmDescript2 = list_rmDescripts{2};
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
rmroiCellSameVox = cell(size(rmroiCell));
for jj = 1:numRois
    for ii = 1:numSubs        
        % get identical voxels for each subject's roi over all ret models
        D = rmroiCell(ii,jj,:);
        % GLU EDIT function: remove voxels from the oppossite hemifield
        rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, cr.defaults.covfig.vfc);
    end
end

% Linearize the data
% Take the difference between 2 rms. 
% Also store the x and y data
for jj = 1:numRois
    % initializing the difference of the centers' thetas
    ldata = []; 

    % intializing the location of the centers
    xdata_rm1 = [];
    ydata_rm1 = []; 
    xdata_rm2 = [];
    ydata_rm2 = [];
    
    % initializing eccentrcity
    ecc_rm1   = [];
    ecc_rm2   = []; 
    
    % initializing angle
    ph_rm1    = [];
    ph_rm2    = []; 
    
    % initializing size
    sm_rm1    = [];
    sm_rm2    = []; 
    
    for ii = 1:numSubs
        rmroi1 = rmroiCellSameVox{ii,jj,1};
        rmroi2 = rmroiCellSameVox{ii,jj,2};
        
        % some subjects don't have 
        if ~isempty(rmroi1) & ~isempty(rmroi2)
            data1 = rmroi1.ph;
            data2 = rmroi2.ph;

            % the difference between centers' thetas.
            % this will determine the color of the line
            % we take absolute value because we are interested in the magnitude
            % of the rotation and not the direction
            fieldDiffOver = abs(data2 - data1);  

            % Note that the difference will range between 0 and 2pi. 
            % We want to constrain values to be between and pi (again not 
            % interested in the direction of the rotation but the magnitude of it)
            % For values greater than pi, subtract it from 2pi
            fieldDiff = ff_polarAngleBetween0AndPi(fieldDiffOver);

            ldata = [ldata fieldDiff];

            % the location of the pRF centers
            xdata_rm1 = [xdata_rm1 rmroi1.x0]; 
            ydata_rm1 = [ydata_rm1 rmroi1.y0]; 

            xdata_rm2 = [xdata_rm2 rmroi2.x0]; 
            ydata_rm2 = [ydata_rm2 rmroi2.y0]; 

            ecc_rm1   = [ecc_rm1 rmroi1.ecc];
            ecc_rm2   = [ecc_rm2 rmroi2.ecc];
            
            ph_rm1    = [ph_rm1 rmroi1.ph];
            ph_rm2    = [ph_rm2 rmroi2.ph];
            
            sm_rm1    = [sm_rm1 rmroi1.sigma1];
            sm_rm2    = [sm_rm2 rmroi2.sigma1];
        end        
    end
    L_data{jj} = ldata; 
    
    X_rm1{jj}   = xdata_rm1;
    Y_rm1{jj}   = ydata_rm1;
    
    X_rm2{jj}   = xdata_rm2;
    Y_rm2{jj}   = ydata_rm2;
    
    Ecc_rm1{jj} = ecc_rm1; 
    Ecc_rm2{jj} = ecc_rm2;
    
    Ph_rm1{jj}  = ph_rm1; 
    Ph_rm2{jj}  = ph_rm2;
    
    Sm_rm1{jj}  = sm_rm1; 
    Sm_rm2{jj}  = sm_rm2;
    
end
% Get a colormap according to the linearized data in L_data
for jj = 1:numRois    
    ldata = L_data{jj};     
    cdata = ff_colormapForValues(ldata, cr.defaults.covfig.vfc.cmapValues, ...
                                        cr.defaults.covfig.vfc.cmapRange);    
    C_data{jj} = cdata;  
end

% colormap for histogram
% cmapValuesHist = colormap('pink');
% cmapValuesHist_tem = colormap('hot');
% cmapValuesHist = cmapValuesHist_tem(2:55, :); 
colormap(zeros(64,3)); % matlab has funky behavior where the size of this influences the size of all future colorbars...
cmapValuesHist = colormap('pink');
close; 

%% PLOT IT
for fov      = [ 1.5]
xx = mrvNewGraphWin('LineRadiality and Scatterplot','wide');
set(xx,'Position',[0.005 0.062 .95 .55 ]);
for jj = 1:numRois  
    % data
    ldata = L_data{jj}; 
    
    X1 = X_rm1{jj};
    Y1 = Y_rm1{jj};
    
    X2 = X_rm2{jj};
    Y2 = Y_rm2{jj};
    
    C = C_data{jj};
    
    ecc_rm1 = Ecc_rm1{jj}; 
    ecc_rm2 = Ecc_rm2{jj};
    
    ph_rm1 = Ph_rm1{jj}; 
    ph_rm2 = Ph_rm2{jj};
    
    sm_rm1 = Sm_rm1{jj}; 
    sm_rm2 = Sm_rm2{jj};
    
    roiName = list_roiNames{jj};
    
  
    % plot on polar map
    %   + black if both word and cb are in the same quadrant (); else gray
    %   + continuous line if word is inside fovea (0-x (x=2deg)), and CB is out; else dashed
    %   + count percentage of continuous black lines over rest:
    %      + bin by size (small sizes very unreliable)
        
    % Filter by size, always, +-inf for no filter
    % Moved this threshold above, now I can't count, make reports in the f
    % sizeMIN  = .5;
    % sizeMAX  = 8;
    
    % sIndw = (sm_rm1 <= sizeMAX) & (sm_rm1 >= sizeMIN); 
    % sIndc = (sm_rm2 <= sizeMAX) & (sm_rm2 >= sizeMIN); 
    % sInd  = (sIndw & sIndc);
    
    % ldata = ldata(sInd); 
    % X1 = X1(sInd);
%     Y1 = Y1(sInd);
%     X2 = X2(sInd);
%     Y2 = Y2(sInd);
%     C = C(sInd);
%     ecc_rm1 = ecc_rm1(sInd); 
%     ecc_rm2 = ecc_rm2(sInd);
%     ph_rm1  = ph_rm1(sInd); 
%     ph_rm2  = ph_rm2(sInd);
%     sm_rm1  = sm_rm1(sInd); 
%     sm_rm2  = sm_rm2(sInd);
%     % histogram(sm_rm1);hold on;histogram(sm_rm2);legend()
    % Obtain angle and ecc again, just in case
    % [PW,ECCW] = cart2pol(X1,Y1);
    % [PC,ECCC] = cart2pol(X2,Y2);
    
   if 0
        fovs   = [0,0.5,1,1.5, 2,3,4,5];
        for ff = 1:length(fovs)
            fov = fovs(ff);
            counter = 0;
            for pp = 1:length(X1)
                pw   = rad2deg(PW(pp));
                pc   = rad2deg(PC(pp));
                qw   = floor(pw/90)+1;
                qc   = floor(pc/90)+1;
                eccw = ECCW(pp);
                eccc = ECCC(pp);

                % if qw==qc && eccw < fov && eccc > fov; counter=counter+1; end
                % if qw==qc && (eccc - eccw) > fov; counter=counter+1; end
                % if qw==qc && (eccc > eccw); counter=counter+1; end
                if (-eccc + eccw) > fov; counter=counter+1; end
                % if eccw < fov && eccc > fov; counter=counter+1; end
                % if qw==qc && ((eccc - eccw) > fov) ; counter=counter+1; end
            end
            if ff==1
                fprintf('\n\n--------------------------------------------------\n');
                % fprintf('SizeMIN: %g, sizeMAX: %g', sizeMIN, sizeMAX);
                % fprintf('  (orig: %g, filtered: %g)', length(sInd), sum(sInd));
                fprintf('\n---------------------------------------------------\n');
            end
            fprintf('Same quadrant and word inside %g deg and cb out, %i out of %i, %0.2f%%\n',...
                    fov,counter,length(X1),counter*100/length(X1))
        end
   end
    testing = false;
    radialityCone=true;
   if 1
        % initialize polar plot
        % Limit plot to visual field circle
        subplot(1,numRois,jj); 
        axis([-cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange ...
              -cr.defaults.covfig.vfc.fieldRange cr.defaults.covfig.vfc.fieldRange])

        % polar plot
        ff_polarPlot(cr.defaults.covfig.vfc); 
        hold on; 

        % colorbar
        % c = colorbar;
        % colormap(cr.defaults.covfig.vfc.cmapValues)
        % set(c, 'Color', [1 1 1])
        % caxis(cr.defaults.covfig.vfc.cmapRange)
        counter  = 0;
        
        maxang15 = [];
        minang15 = [];
        maxang5  = [];
        minang5  = [];
        for pp = 1:length(X1)
            % pp = 555;
            % lineColor = C(pp,:);
            % pw   = rad2deg(PW(pp));
            % pc   = rad2deg(PC(pp));
            % qw   = floor(pw/90)+1;
            % qc   = floor(pc/90)+1;
            % eccw = ECCW(pp);
            % eccc = ECCC(pp);
            
            pw   = ph_rm1(pp);
            pc   = ph_rm2(pp);
            eccw = ecc_rm1(pp);
            eccc = ecc_rm2(pp);
            
            

            % if qw==qc ; lineColor = [1 1 1];
            % else; lineColor = [.5 .5 .5]; end

            % if eccw < 5 && eccc > 5 ; lineStyle = '-';
            % else; lineStyle = ':'; end

            % if qw==qc && eccw < 5 && eccc > 5; counter=counter+1; end
            if (abs(eccc - eccw)  > fov)
                if  ((eccc - eccw) > 0) % qw==qc &&
                    counter=counter+1; 
                    lineColor = [.1 .1 .1];
                    lineStyle = ['-'];
                else
                    lineColor = [.5 .5 .5];
                    lineStyle = [':'];
                end
                
                if radialityCone
                    % WORDS
                    x1 = X1(pp);
                    y1 = Y1(pp);
                    % CB 
                    x2 = X2(pp);
                    y2 = Y2(pp);
                    % We know the polar angle and the ecc, rotate & translate
                    % Rotate angle
                    pcn    = 0;
                    pwn    = pw - pc;
                    % Convert to cartesian
                    [x1n, y1n] = pol2cart(pwn, eccw);
                    [x2n, y2n] = pol2cart(pcn, eccc);
                    % Do the translation only in the x axis
                    if eccc > 5
                        x2nt  = 15;
                        x2dif = x2nt - x2n;
                        x1nt  = x1n + x2dif;
                        % Plot with the same y
                        line([x1nt x2nt], [y1n, y2n], ...
                            'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',1); 
                    else
                         x2nt  = 5;
                         x2dif = x2nt - x2n;
                         x1nt  = x1n + x2dif;
%                         % Plot with the same y
%                         line([x1nt x2nt], [y1n, y2n], ...
%                             'Color', lineColor, 'LineStyle', lineStyle, 'LineWidth',1); 
                    end
                    
                     % Calculate angles
                     angle = atand(y1n/(x2nt-x1nt));
                    if  ((eccc - eccw) > fov)
                        if eccc > 5
                            if angle > 0
                                maxang15 = [maxang15 angle];
                            end
                            if angle < 0
                                minang15 = [minang15 angle];
                            end
                        else
                            if angle > 0
                                maxang5 = [maxang5 angle];
                            end
                            if angle < 0
                                minang5 = [minang5 angle];
                            end
                        end
                    end
                    
                    if testing
                        line([x1 x2], [y1, y2], ...
                            'Color', 'r', 'LineStyle', lineStyle, ...
                             'LineWidth',1); 
                         % Assertions that the calculations are right
                        d      = sqrt((x2-x1)^2 +  (y2-y1)^2)
                        dn     = sqrt((x2nt-x1nt)^2 +  (y2n-y1n)^2)
                    end
                    
                else
                    plot([X1(pp) X2(pp)], [Y1(pp), Y2(pp)], ...
                        'Color', lineColor, 'LineStyle', lineStyle, ...
                         'LineWidth',1);       
                end
            end
        end

        % titleName = {
            % ['pRF Center radiality in ' roiName]
            % ['white(' sprintf('%0.2f%% ',counter*100/length(X1)) ...
            %  ')=same quadrant, ecc: (' rmDescript2 ' - ' rmDescript1 ') > ' num2str(fov)]
            % ['sMIN:' num2str(sizeMIN) ',sMAX:' num2str(sizeMAX) ...
            % ' (pRFs ' num2str(sum(sInd)) ' of '  num2str(length(sInd)) ')']
        %     };
        titleName = {
            [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_')]; 
            % ['slope: ' num2str(meanSlope)];
            % ['ci: ' num2str(ci')];
            % [num2str(numVoxels) ' voxels']
            };
        ttext = strrep(titleName{1},'\_left','');
        
        % title(titleName, 'fontweight', 'bold', 'color', [.1 .1 .1], 'fontsize', 14);
        text(-9, 11, ttext,'Color','k','FontWeight','Bold','FontSize',16)
        
        if ~isempty(maxang15)
            text(-9, 3,sprintf('+%1.1f°',max(maxang15)),'Color','k','FontSize',10)
            text(-9, 1,sprintf('+%1.1f°',median(maxang15)),'Color','r','FontSize',10)
            line([15+5*cosd(180-median(maxang15)),15],[5*sind(median(maxang15)),0],'Color','r')
        end
        if ~isempty(minang15)
            text(-9, -1,sprintf('%1.1f°',median(minang15)),'Color','r','FontSize',10)
            text(-9, -3,sprintf('%1.1f°',min(minang15)),'Color','k','FontSize',10)
            line([15+5*cosd(180-median(minang15)),15],[5*sind(median(minang15)),0],'Color','r')
        end
        
%         if ~isempty(maxang5)
%             text(-9, 3,sprintf('+%1.1f°',max(maxang5)),'Color','k','FontSize',10)
%             text(-9, 1,sprintf('+%1.1f°',median(maxang5)),'Color','r','FontSize',10)
%             line([5+5*cosd(180-median(maxang5)),5],[5*sind(median(maxang5)),0],'Color','r')
%         end
%         if ~isempty(minang5)
%             text(-9, -1,sprintf('%1.1f°',median(minang5)),'Color','r','FontSize',10)
%             text(-9, -3,sprintf('%1.1f°',min(minang5)),'Color','k','FontSize',10)
%             line([5+5*cosd(180-median(minang5)),5],[5*sind(median(minang5)),0],'Color','r')
%         end
        xlim([-15, 18])
        % titlefile = strrep(titleName{1},' ','_');
        % saveas(gcf, fullfile(crRootPath,'local','png',[titlefile '.png']), 'png') 
%         subplot(2,numRois,numRois+jj); 
%         polarhistogram(deg2rad([maxang5 minang5]),100, 'Normalization','probability',...
%             'EdgeAlpha',1, 'EdgeColor','k','FaceAlpha',1, 'FaceColor','k'); 
%         hold on;
%         polarhistogram(deg2rad([maxang15 minang15]),100, 'Normalization','probability',...
%             'EdgeAlpha',1, 'EdgeColor','b','FaceAlpha',.4, 'FaceColor','b'); 
%         legend({'< 5°','> 5°'})
        

   end
end
saveas(gcf, fullfile(crRootPath,'local','png',['Radiality_Plots_2x' num2str(fov) '.png']), 'png')    
saveas(gcf, fullfile(crRootPath,'local','svg',['Radiality_Plots_2x' num2str(fov) '.svg']), 'svg')    
close all
end

%% SCATTERPLOTS


% whether looking at a subject by subject basis
subIndividually = false; 

% field to plot. Ex:  
% 'co': variance explained 
% 'ecc': eccentricity 
% 'sigma': effective size 
% 'sigma1': sigma major
% 'numvoxels' for number of voxels in roi
% fieldToPlotDescript is for axis labels and plot titles
%     'sigma'       : effective sigma
%     'sigma1'      : sigma major
%     'ecc'         : eccentricity
%     'co'          : variance explained 
%     'exponent'    : exponent
%     'betaScale'   : how much to scale the predicted tseries by
%     'meanMax'     : mean of the top 8 values
%     'meanPeaks'   : mean of the outputs of matlab's meanPeaks
list_fieldNames  = {    
    'co'
    'ecc'
%     'sigma'
%     'ph'
    }; 

list_fieldDescripts = {
    'variance explained'
    'eccentricity'
%     'sigma'
%     'polar angle'
    }; 

% which plots do we want? lots we can make ...
plot_fit = false; % plotting the across-subject bootstrapped line w/ CIs

% transparency of the plots
alphaValue = 0.4; 



% location of the colorbar
% default: 'eastoutside'
% 'southoutside': 
cbarLocation = 'eastoutside';

% end modification section

numSubs = length(list_subInds);
numRois = length(list_roiNames);
numRms = length(list_rmNames);

% number of fields
numFields = length(list_fieldNames);

% rm descriptions
rm1Descript = list_rmDescripts{1}; 
rm2Descript = list_rmDescripts{2}; 

% initialize structs / matrices for mixed effects
subjectLines = cell(numSubs, numRois, numFields); % because pairwise

% initialize struct for calculating the percentage of voxels above the
% identity line
percentAbovePooled = zeros(numRois, numFields);
percentAboveSubs   = zeros(numSubs, numRois, numFields);
A = cell(numFields*numRois, 5);

% get the cell of rms so that we can threshold
% rmroiCell = ff_rmroiCell(cr, list_subInds, list_roiNames, list_dtNames, list_rmNames, ...
%     'list_path', list_path);

% Threshold and get identical voxels for each subject
% In comparing ret models, the collection of voxels may not be the same
% because of the thresholding. In this cell we redefine the rmroi
% rmroiCellSameVox = cell(size(rmroiCell));

% for jj = 1:numRois
%     for ii = 1:numSubs        
        % get identical voxels for each subject's roi over all ret models
%         D = rmroiCell(ii,jj,:);
%         rmroiCellSameVox(ii,jj,:) = ff_rmroiGetSameVoxels(D, vfc);        
%     end
% end

% close all;


    % subplot(2,numRois,numRois+jj); 

    
for ff = 1:numFields

    % field-specific properties
    fieldName = list_fieldNames{ff};
    fieldNameDescript = list_fieldDescripts{ff}; 

    if strcmp(fieldName, 'sigma1') 
        maxValue = cr.defaults.covfig.vfc.sigmaMajthresh(2);
    elseif strcmp(fieldName, 'sigma')          
        maxValue = cr.defaults.covfig.vfc.sigmaEffthresh(2);
    elseif strcmp(fieldName, 'ecc')
        maxValue = cr.defaults.covfig.vfc.eccthresh(2);
        fov = 1.5; % width of the band
        nrows = 2; ncols = 3;
        position = [0.005 0.062 .9 .7 ];
    elseif strcmp(fieldName, 'co')
        maxValue = 1; 
        fov = 0.2; % width of the band
        nrows = 1; ncols = 6;
        position = [0.005 0.062 .9 .5 ];
    elseif strcmp(fieldName, 'exponent')
        maxValue = 2; 
    elseif strcmp(fieldName, 'meanMax')
        maxValue = 20; 
    elseif strcmp(fieldName, 'meanPeaks')
        maxValue = 10; 
    elseif strcmp(fieldName, 'betaScale')
        maxValue = 5; 
    elseif strcmp(fieldName, 'x0') || strcmp(fieldName, 'y0')
        maxValue = cr.defaults.covfig.vfc.fieldRange;
    elseif strcmp(fieldName, 'ph')
        maxValue = 2*pi; 
    else
        error('Define the maxValue so we can normalize and fit the beta distribution.');
    end
    
    xx = mrvNewGraphWin('Scatterplots');
    set(xx,'Position',position);
    for jj = 1:numRois
        roiName = list_roiNames{jj};
        subplot(nrows,ncols,jj);
    

        
        axisLims = [0 maxValue]; 
        BarData1 = [];
        BarData2 = [];

        for ii = 1:numSubs

            subInd = list_subInds(ii);
            
            % rmRois for different ret models
            rmroi1 = rmroiCellSameVox{ii,jj,1}; 
            rmroi2 = rmroiCellSameVox{ii,jj,2};

             if ~isempty(rmroi1)
                 
                % get the data
                x1 = eval(['rmroi1.' fieldName]);
                x2 = eval(['rmroi2.' fieldName]);

                % shift so that the smallest value is 0
                if strcmp(fieldName, 'ph')
                    x1 = x1 + pi; 
                    x2 = x2 + pi;
                end
                
                %% For the mixed effects
                % fit a line for each subject
                % p = polyfit(x1, x2, 1);              
                % subjectLines{ii,jj,ff} = p;                 
                % [B, BINT] = regress(x2', x1'); 
                % b.slope     = B; 
                % b.ci95      = BINT; 
                
                % p = polyfit(x1, x2, 1);
                % b.pintercept = p(1);
                % b.pslope = p(2); 
                
                % subjectLines{ii,jj,ff} = b; 
                
                %% the percentage of voxels above the identityLine
                perAbove = sum(x2 > x1) / length(x2);
                percentAboveSubs(ii,jj,ff) = perAbove; 
                
                %% concatenate
                BarData1 = [BarData1, x1];
                BarData2 = [BarData2, x2];
           
            end
        end % end loop over subjects
        
        %% mixed effects: fit a line to individual subjects
        % slopes = nan(1, numSubs);   
        % slopesp = nan(1,numSubs);
        % interceptsp = nan(1,numSubs);
        % percents = percentAboveSubs(:,jj,ff);
        
        % for ii = 1:numSubs
        %     b = subjectLines{ii,jj,ff}; 
        %     if ~isempty(b)
        %         slopes(ii) = b.slope; 
        %     end
        % end

        % the calculating. nan will cause bootci to error
        % slopes(isnan(slopes)) = []; 
        % percents(isnan(percents)) = [];
        
        % table things. 
        % (1)roiName (2)fieldName (3) ciLow (4)ciHigh (5)mean
        % tind = (jj-1)*numFields + ff; 
        
        
        % if numSubs > 1
            % numbs = 1000; 
            % [ci, bootstat] = bootci(numbs, @mean, slopes);
            % meanSlope = mean(bootstat); 
            
            % % table things for percent above identityLine
            % [ciPer, bootstatPer] = bootci(numbs, @mean, percents);
            % A{tind,1} = ff_stringRemove(roiName, '_rl');
            % A{tind,2} = fieldName;
            % A{tind,3} = ciPer(1);
            % A{tind,4} = ciPer(2);
            % A{tind,5} = mean(bootstatPer);
            
        % else
        %     meanSlope = nan; 
        %     ci = nan; 
        %   end
                
        %% calculations related to the percentage above the identityLine
        % percentabove = sum(BarData2 > BarData1) / length(BarData1);
        % percentAbovePooled(jj,ff) = percentabove;
        
        % properties related to both types of scatter plots
        % coloring by number of voxels or percentage of voxels
        npoints = 100; 
        
        % 3d histogram heat map -- absolute number of voxels
        ff_histogramHeat(BarData1, BarData2, maxValue, maxValue, 50,cmapValuesHist,fov,roiName);
        numVoxels = length(BarData1); 
        % axes and title
        switch fieldName
            case {'ecc'}
                if jj==1 || jj==4
                    ylabel(['pRF eccentricity for ' rm2Descript ' (deg)'])
                end
                if jj>3
                    xlabel(['pRF eccentricity for ' rm1Descript ' (deg)'])
                end
            case {'co'}
                if jj==1 % || jj==4
                    ylabel(['Variance explained for ' rm2Descript ' (%)'])
                end
                % if jj>3
                    xlabel(['Variance explained for ' rm1Descript ' (%)'])
                % end
            otherwise
                if jj==1 || jj==4
                    ylabel(['' rm2Descript ''])
                end
                if jj>3
                    xlabel(['' rm1Descript ''])
                end
        end
        
    end % loop over rois

    % titleName = {
    %     [strrep(ff_stringRemove(roiName, 'WangAtlas_'),'_','\_') '.' fieldName];
    % ['slope: ' num2str(meanSlope)];
    % ['ci: ' num2str(ci')];
    % [num2str(numVoxels) ' voxels']
    %      };
    % title(titleName, 'FontWeight', 'Bold');
    saveas(gcf, fullfile(crRootPath,'local','png',[titleName{1} '_' fieldName '_band-2x' num2str(fov) '.png']), 'png')
    saveas(gcf, fullfile(crRootPath,'local','svg',[titleName{1} '_' fieldName '_band-2x' num2str(fov) '.svg']), 'svg')
end % loop over fields

% percent above identityLine, pooled over subjects ... print out the 
% percentAbovePoooled

% Percent above identity line. Bootstrapped across subjects (mixed effects)
% {   ['Percent of voxels above identity line. ']
%     ['Bootstrapped across subjects']
%     [rm2Descript ' vs. ' rm1Descript]
% }
% T = cell2table(A, 'VariableNames', {'roiName', 'fieldName', 'ciLow', 'ciHigh', 'MeanPercent'});

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






=======
%}









