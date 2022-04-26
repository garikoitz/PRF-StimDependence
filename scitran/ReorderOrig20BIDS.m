% This scripts will dow


%% 1.- 

% Read the bookkeeping file from Rosemary too
cr         = struct();
cr.codeDir = crRP;

% WHERE THE NEW DATA IS
cr.dirs.BASE     = '/black/localhome/glerma/TESTDATA/PRF-StimDependence';
cr.dirs.DATA     = fullfile(cr.dirs.BASE,'DATA');
cr.dirs.ANALYSIS = fullfile(cr.dirs.BASE,'ANALYSIS');
cr.dirs.BIDS     = fullfile(cr.dirs.BASE,'BIDS');
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
% addpath(genpath(fullfile(cr.dirs.ANALYSIS,'matlabfiles')));

% Rosemary relied on this file that contains most of the subjects and other
% lists. Make it work with relative paths and store it in each project repository
% This file was used as well for: 
% - copying files to a new location
% - editing mrSession to reflect the file changes
cr.bk = bookKeeping(cr);

%% 2.- Create table with files within each session
% Find the subjects we are interested in FW
% list_subInds        = [31:36 38:44];  % Hebrew
list_subInds        = [1:20];  % CNI 20
TR = 2;

ses = 'T01';


ii = 0;
for subind = list_subInds
    ii = ii + 1;
    subname = cr.bk.list_sub{subind};
    [~,anatName]=fileparts(cr.bk.list_anatomy{subind});
    fprintf('\nSubDetails:\nInd:%i, StrInd:%s, subname:%s, Name:%s, anatName:%s\n',...
        subind,cr.bk.list_subNumberString{subind},subname,...
        cr.bk.list_names{subind},anatName)
    bidsname = sprintf('S%03i', subind);
    % Read the init file to know where the files are: 
    a = fileread(fullfile(cr.bk.list_sessionRet{subind},'pp_mrInit.m'));
    b = strrep(a,'mrInit(params)','%mrInit(params)');
    c = strrep(b,'clear all; close all; clc;','%clear all; close all; clc;');
    d = strrep(c,'% params.sliceOrder','params.sliceOrder');
    e = strrep(d,'check_exist','% check_exist');
    % Run it to obtain all the paths we want
    eval(e)
    % Check that params.sliceOrder exists
    if ~isfield(params,'sliceOrder');warning('For %i sliceOrder does not exist',subind);end
  
    % BIDSIFY THE FILES
    % anat first
    src  = path_anatomy;
    dstd = fullfile(cr.dirs.BIDS,['sub-' bidsname],['ses-' ses],'anat');
    dstf = fullfile(dstd, ['sub-' bidsname '_ses-' ses '_T1w.nii.gz']);
    [~,d,m] = myCopy(src,dstf,subind,false);
    
    % Now the functional files
    dstfs = {}; ii = 0;
    for np =1:length(path_functionals)
        [p,f,e] = fileparts(path_functionals{np};km);
        if ~contains(p,'Localizer_')
            ii = ii +1;
            src  = fullfile(p,'func.nii.gz');
            % Obtain the task
            pp = split(p,"/");
            task = char(strrep(pp{end},"Ret_",""));
            run  = task(end);
            task = task(1:end-1);
            dstd = fullfile(cr.dirs.BIDS,['sub-' bidsname],['ses-' ses],'func');
            dstf = fullfile(dstd, ['sub-' bidsname '_ses-' ses '_task-' task '_run-' run '_bold']);
            dstfn = [dstf '.nii.gz']; 
            [s,d,m] = myCopy(src,dstfn,subind,false);
            % Now write the json file with the basic info
            st = TR * (params.sliceOrder - 1) / length(params.sliceOrder);
            js = struct('SliceTiming',st);
            js.RepetitionTime = 2;
            js.Modality       = 'MR';
            js.MagneticFieldStrength = 3;
            jsonwrite([dstf '.json'], jsonencode(js));
            dstfs{ii} = dstf;
        end
    end
    % Now the stimulus files
    % First copy the files as they are right now:
    % IMAGE file
    imFile_Knk          = fullfile(fileparts(p),'Stimuli/images_knk_fliplr.mat');              % Words and FalseFont
    [~,fn] = fileparts(imFile_Knk);
    dstfiknk  = fullfile(cr.dirs.BIDS,'sourcedata','stimuli',[fn '.mat']);
    [s,d,m] = myCopy(imFile_Knk,dstfiknk,subind,false);
    
    imFile_Checkers     = fullfile(fileparts(p),'Stimuli/images_8barswithblank_fliplr.mat');   % Checkers
    [~,fn] = fileparts(imFile_Checkers);
    dstfichck  = fullfile(cr.dirs.BIDS,'sourcedata','stimuli',[fn '.mat']);
    [s,d,m] = myCopy(imFile_Checkers,dstfichck,subind,false);
    
    % PARAMS file
    paramsFile_Knk      = fullfile(fileparts(p),'Stimuli/params_knkfull_multibar_blank.mat');  % Words and FalseFont
    % Add the loadMatrix field to the params files
    tt = matfile(paramsFile_Knk,'writable',true);
    pms = tt.params; 
    pms.loadMatrix = imFile_Knk;
    tt.params = pms;
    clear tt; clear pms;
    [~,fn] = fileparts(paramsFile_Knk);
    dstfpknk  = fullfile(cr.dirs.BIDS,'sourcedata','vistadisplog',['sub-' bidsname],['ses-' ses],[fn '.mat']);
    [s,d,m] = myCopy(paramsFile_Knk,dstfpknk,subind,false);
   
    paramsFile_Checkers = fullfile(fileparts(p),'Stimuli/params_checkers.mat');  % Checkers
    % Add the loadMatrix field to the params files
    tt = matfile(paramsFile_Checkers,'writable',true);
    pms = tt.params; 
    pms.loadMatrix = imFile_Checkers;
    tt.params = pms;
    clear tt; clear pms;
    [~,fn] = fileparts(paramsFile_Checkers);
    dstfpchck  = fullfile(cr.dirs.BIDS,'sourcedata','vistadisplog',['sub-' bidsname],['ses-' ses],[fn '.mat']);
    [s,d,m] = myCopy(paramsFile_Checkers,dstfpchck,subind,false);
    
    
    
    % Create one symlink per each task and run
    cd(fileparts(dstfpknk))
    for nd = 1: length(dstfs)
        [~,dstf] = fileparts(dstfs{nd});
        dstf = strrep(dstf,"1_bold","01_params.mat");
        dstf = strrep(dstf,"2_bold","02_params.mat");
        dstf = strrep(dstf,"3_bold","03_params.mat");
        dstf = strrep(dstf,"4_bold","04_params.mat");
        % dstf = fullfile(fileparts(dstfpknk),dstf);
        ppp  = split(dstf,"_");
        task = split(ppp(3),"-"); task = task{2};
        switch task
            case {'English','Hebrew'}
                if ~isfile(char(dstf))
                    [~,f,e] = fileparts(dstfpknk);
                    ss = system(['ln -s ./' f e ' ./' char(dstf)]);
                end
            case {'Checkers'}
                 if ~isfile(char(dstf))
                     [~,f,e] = fileparts(dstfpchck);
                    system(['ln -s ./' f e ' ./' char(dstf)]);
                 end
            otherwise
                warning('Task %s not recognized',task)
        end    
    end
    
    

end

% Create the required BIDS files
% I created them manually, they are in the base of BIDS_HEB

% Fix the jason files, they have a weird \" instead of just \






