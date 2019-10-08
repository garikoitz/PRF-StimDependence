function cr = cr_mrInit(cr,opt,varargin)
%function for mrInit: mrVista initialization
% specifies paths and checks that they exist
% clips frames
% runs motion correction
%
% rl, summer 2014 
% glu, summer 2019: took rl's script, and make into a function
% objectives: 
%   - reproduce rl's results
%   - make it generic so that we can use it with other projects. 

% clear all; close all; clc;

%% Read parameters (this is usually the modify section of RL's scripts)

% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('cr'  , @isstruct);
p.addRequired('opt' , @isstruct);

% Parse. Assign result inside each case
p.parse(cr, opt, varargin{:});
% Read here only the generic ones
% opt = p.Results.opt;


% Add other lists coming from bookKeeping
% The following parameters have been set in bk = bookKeeping():
% session list. see bookKeeping
opt.list_path            = cr.bk.list_sessionRet; 
% list of checker scan number, corresponding to session list
opt.list_numCheckers     = cr.bk.list_scanNum_Checkers_sessionRet; 
% list of knk scan number, corresponding to session list
opt.list_numKnk          = cr.bk.list_scanNum_Knk_sessionRet;  


% create params structure ...
params = mrInitDefaultParams; 

%% Do the initialization
%% loop over all subjects
for ii = opt.list_subInds
% specify session path. usually this script is saved right there 
    cr.dirVista = opt.list_path{ii};
    chdir(cr.dirVista);
    
    % specify inplane
    % TODO: where is this generated? should come from dicoms
    cr.path_inplane = fullfile(cr.dirVista, 'prescribeInplane1','inplane_xform.nii.gz'); 
    
    % specify 3DAnatomy file
    cr.path_anatomy = fullfile(cr.bk.list_anatomy{ii}, 't1.nii.gz');

    % note for each of the functional scans
    params.annotations  = {
        'Localizer_Hebrew1'
        'Localizer_Hebrew2'
        'Ret_Checkers1'
        'Ret_English1'
        'Ret_Hebrew1'
        'Ret_Hebrew2'
        };
    
    % specify the functional files
    % TODO: pass this as parameter as well
    cr.path_functionals = {};
    for na = 1:length(params.annotations)
        cr.path_functionals = [cr.path_functionals; ...
                               {fullfile(cr.dirVista, params.annotations{na}, ...
                                            'func_xform.nii.gz')}];
    end
    % subject name
    params.subject      = lower(strrep(cr.bk.list_names{ii},' ',''));

    % session code
    % params.sessionCode  = 'heb_pilot09_avbe'; 
    parts = split(cr.dirVista,filesep);
    params.sessionCode  = [parts{end-1} '_' params.subject];
                       
    % description
    params.description  = ['ret: ' parts{end}]; 
    
    % specify frames to keep. nScans x 2 matrix describing which frames to keep from 
    % each scan's time series. The first column specifies the number to skip
    % the second column specifies the number to keep. A flag of -1 in the 2nd
    % column indicates to keep all after the skip. Leaving everything empty
    % will cause all frames to be kept (the default)
    
    % TODO: this needs to come from bookKeeping
    params.keepFrames = [
        [4, 93];
        [4, 93];
        [10, 144];
        [10, 144];
        [10, 144];
        [10, 144];
        ];
    % specify parfiles (for localizer) {1 x nScans}
    % paths can be absolute or relative to Stimuli/parfile
    % params.parfile = {
    %     [];
    %     [];
    %     [];
    %     };
    
    
    %%  automated tests to check everything is in place
    
    % check that all specified files exist
    check_exist(cr.dirVista);
    check_exist(cr.path_inplane);
    check_exist(cr.path_anatomy);
    check_exist(cr.path_functionals);
    
  
    % ... and insert the required parameters
    params.sessionDir   = cr.dirVista;
    params.inplane      = cr.path_inplane;
    params.vAnatomy     = cr.path_anatomy;
    params.functionals  = cr.path_functionals;
    
    % check that length of annotions matches number of functionals
    if (length(params.annotations) ~= length(cr.path_functionals))
        error('Check that annotionals align with functionals');
    end
    
    % check that length of keep frames matches number of functionals
    if (size(params.keepFrames,1) ~= length(cr.path_functionals))
        error('Check that annotionals align with functionals');
    end
    
    % motion compensation.
    % 4 means between and within, with within being first
    params.motionComp = 4;
    
    % slice timing correction/
    % params.sliceTimingCorrection = 1;
    % params.sliceOrder = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36];
    
    %% Go!
    mrInit(params); 
    % Return params in cr
    cr.mrInitParams = params;
    
end






