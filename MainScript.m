%% STEPS
% This repository is based on the code created by Rosemary Le, part of coverageReading
% repository. For the Stimulus Dependence paper (2019) we tried to make it into
% a reproducible and reusable process. The thing is that we need to separate
% this specific project from others, we need to separate data from code. 
close all; clear all;
mrvCleanWorkspace;
cr         = struct();
cr.codeDir = crRootPath;
cr.baseDir = '/Users/glerma/Documents/STANFORD_PROJECTS';
cr.dataDir = fullfile(cr.baseDir, 'PRF-StimDependence');
% add to path the required matlab files with info to run the project
addpath(genpath(fullfile(cr.dataDir, 'matlabFiles')));

% Rosemary relied on this file that contains most of the subjects and other
% lists. Make it work with relative paths and store it in each project
% repository
cr.bk = bookKeeping(cr.dataDir);
% subjects we want to do this for
opt.list_subInds        = [31]; 
chdir(cr.bk.list_sessionRet{opt.list_subInds});

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

% apply canonical x form -- for every nifti
niftiApplyCannonicalXform

% acpc align the anatomical 
mrAnatAverageAcpcNifti

% run freesurfer
eval(['! recon-all -i ' pathT1 ' -subjid ' dirNameFreesurfer ' -all'])

% ribbon from freesurfer into class file -- t1_class.nii.gz
fs_ribbon2itk(inputRibbonFile, outputClassNii, [], pathT1, [])

% TODO: see if we can substitute the previous steps using fmriprep


%% Initialization 
% Initialize. This will create the mrSESSION in the root folder. 
cr_mrInit(cr, opt);


%% Allignment and dataType creation
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

%% Stimuli preparation
% make a Stimuli folder in the same place as the mrSESSION.mat
% for localizer GLM analyses, make Stimuli/Parfiles
% 2 things go into Stimuli
% params file -- mrVista writes this to desktop
% image matrix -- (part of the params file)



%% Analysis of the data
% The code here will go to a Docker container. Make the figures in the paper reproducible
% The code below will come as a combination from pmVistasoft.m that I did
% adapting a script from Jon Winawer, and s_prfRun.m Rosemary Le's script. 

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


%% Generate the individual figure
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
opt.list_path        = cr.bk.list_sessionRet; 
% ROIs
opt.list_roiNames    = {'lVOTRC'};
% dt and rm names
opt.list_dtNames     = {'Words_Hebrew','Words_English'};
opt.list_rmNames     = {'retModel-Words_Hebrew-css.mat'
                        'retModel-Words_English-css.mat'};
opt.list_rmDescripts = {'Words_Hebrew', 'Words_English'};
cr = figFunction_coverage_individual(cr, opt);
%}

