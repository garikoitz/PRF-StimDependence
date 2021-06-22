function cr = cr_prfRun(cr, subind, varargin)
% Run PRF models
% This was a script from Rosemary. I am going to create a function out of it. 
% It shuold be very similar to the pmVistasoft.m script. We will test both later
% on. 

% GLU: many editions to be able to read all the specifics for each rmName
% and sub, created the folder params to store all the specifics per subject


% 
% 
% Rosemary's notes
% assumes that the datatypes with the averaged tseries are already created
% and xformed
%
% rl  08/2014
% glu 08/2019


% clear all; close all; clc; 
% bookKeeping; 

%% Read parameters (this is usually the modify section of RL's scripts)


% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('cr'      , @isstruct);
p.addRequired('subind'  , @isnumeric);

% Parse. Assign result inside each case
p.parse(cr, subind, varargin{:});
% Read here only the generic ones
% opt = p.Results.opt;


% Individual data
subname         = cr.bk.list_sub{subind};
opt             = cr.subj.(subname).params.prf;
opt.p           = cr.defaults.prfrun.p;
opt.params      = cr.defaults.prfrun.params;

% The following parameters have been set in bk = bookKeeping():
% session list. see bookKeeping


%opt.list_path   = cr.bk.list_sessionRet; 
opt.list_path   = cr.bk.list_sessionSizeRet; 


opt.list_rois   = opt.rois;




% We will obtain this directly in each case, no names associated
% list of checker scan number, corresponding to session list
% opt.list_numCheckers     = cr.bk.list_scanNum_Checkers_sessionRet; 
% list of knk scan number, corresponding to session list
% opt.list_numKnk          = cr.bk.list_scanNum_Knk_sessionRet;  

%% 
    
    % directory with ret vista session. move here
    dirVista = opt.list_path{subind};
    chdir(dirVista);
    
    % open the session
    vw = initHiddenGray;

    % need some global variables later
    load mrSESSION; 

    % main anatomy path
    dirAnatomy = cr.bk.list_anatomy{subind};
    
    % ret parameters based on the subject
    % scan number with checkers and knk, for clip frame information
    % opt.p.scanNum_Knk                  = opt.list_numKnk(subind);
    % opt.p.scanNum_Checkers             = opt.list_numCheckers(subind);
    
    
    %% loop over the datatypes
    for kk = 1:length(opt.list_rmName)

        % set current dataTYPE 
        rmName = opt.list_rmName{kk};
        vw = viewSet(vw, 'curdt', rmName); 

        % get the dataType struct
        dtstruct = viewGet(vw, 'dtstruct'); 

        % get the data type number for later
        dataNum = viewGet(vw, 'curdt'); 

        % some variables depend on whether checkers or knk was run  
        % Left below existing Rosemary code, created function
        % addpath(fileparts(cr.bk.list_prfParams{subind}))
        prf = prfRun_params_rmName(cr, subind, rmName);
        % rmpath(fileparts(cr.bk.list_prfParams{subind}))
        
        opt.params.paramsFile = prf.p.paramsFile;
        opt.params.imFile     = prf.p.imFile;
        opt.p.scanNum         = prf.p.scanNum;
        %{
        if length(rmName) > 7 && strcmp(rmName(1:8), 'Checkers')
            opt.params.paramsFile   = opt.p.paramsFile_Checkers; 
            opt.params.imFile       = opt.p.imFile_Checkers; 
            opt.p.scanNum           = opt.p.scanNum_Checkers; 
        else
            opt.params.paramsFile   = opt.p.paramsFile_Knk; 
            opt.params.imFile       = opt.p.imFile_Knk; 
            opt.p.scanNum           = opt.p.scanNum_Knk; 
        end
        %}


        %% getting parameter values for prf model fit ----------------------
        
        
        switch subind
            case{13,17}
                % This is failing, trying next option

                % Well, subjects 13 and 17 (mw and tl), it seems that this works
                % {
                opt.params.nFrames          = viewGet(vw, 'nFrames');  

                opt.params.framePeriod      = viewGet(vw, 'framePeriod');   
                tem.totalFrames             = mrSESSION.functionals(opt.p.scanNum).totalFrames;  
                opt.params.prescanDuration  = (tem.totalFrames - opt.params.nFrames)*opt.params.framePeriod; 
                %}
            otherwise
        
                % {
                barScans = prf.barScans;


                % In the main files Rosemary ws using, she had the code above
                % NOTE: nFrames and framePeriod are generic
                % I made basScans to be the same as scanNum, so I guess hat the
                % code now is different, plot them
                opt.params.nFrames         = mrSESSION.functionals(barScans(1)).nFrames; 
                opt.params.framePeriod     = mrSESSION.functionals(barScans(1)).framePeriod; 
                tem.totalFrames            = mrSESSION.functionals(barScans(1)).totalFrames;  
                opt.params.prescanDuration = (tem.totalFrames - opt.params.nFrames)*opt.params.framePeriod; 
                %}
        end
        
        
        
        
        
        
        
        
        
        
        
        % store it
        dataTYPES(dataNum).retinotopyModelParams = opt.params;

        % save it
        saveSession; 

        %% Put the rm params into the view structure

        vw = rmLoadParameters(vw);  
        % the function rmLoadParameters used to call both rmDefineParameters
        % and rmMakeStimulus. If we do it here so that we can give it arguments
        % outside of the default (eg previously, sigma major and minor would be 
        % identical despite having prfModel = {'one oval gaussian'} when
        % specifying it as an argument in vw = rmMain(vw, ...)

%         scan/stim and analysis parameters
%         params = rmDefineParameters(vw, varargin)
%         params = rmDefineParameters(vw, 'model', prfModel);
% 
%         make stimulus and add it to the parameters
%          params = rmMakeStimulus(params, keepAllPoints)
%         params = rmMakeStimulus(params);

        % store params in view struct
        vw  = viewSet(vw,'rmParams',opt.params);

        % check it 
        % rmStimulusMatrix(viewGet(vw, 'rmparams'), [],[],2,false);
        
        %% RUN THE PRF!
        %% If we've defined a list of rois to run over, do that loop here
        if ~isempty(opt.list_rois)
            
            for jj = 1:length(opt.list_rois)
                
                % load the current roi
                roiName = opt.list_rois{jj};
                % sigh... anatomy was in two different places, with
                % different names even... ROIs where in one place and fs in
                % another. I copied all fs directory to our project, and
                % then used a matlab script to copy only the ROIs folder
                % and other loose files at the same level where the folders
                % where, just in case. Well, Matlan copied all the files
                % inside all the subfolders, so now all ROIs are loose
                % under the anatomy/subject folder. I comment the original
                % line and edit it removing the ROIs folder
                % roiPath = fullfile(dirAnatomy, 'ROIs', roiName);
                roiPath = fullfile(dirAnatomy, roiName);
                vw = loadROI(vw, roiPath, [], [], 1, 0);
                
                % name the ret model
                outFileName = ['retModel-' rmName '-' opt.prfModel{1} '-' roiName];
                
                % run the model!
                vw = rmMain(vw, [], opt.wSearch, 'model', opt.prfModel, ...
                            'matFileName', outFileName);
                
            end
            
        else
            
            % name the ret model - whole brain
            outFileName = ['retModel-' rmName '-' opt.prfModel{1}];
            
            % no need to load rois, just run it!
            vw = rmMain(vw, [], opt.wSearch, 'model', opt.prfModel, ...
                        'matFileName', outFileName);

        end
           
    end



end
