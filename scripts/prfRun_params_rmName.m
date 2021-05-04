function prf = prfRun_params_rmName(cr, subind, rmName)
% ret parameters based on the subject


% radius of circle retinotopy in visual angle degrees
prf.p.stimSize    = cr.bk.list_stimSize(subind); 

switch rmName
    case {'Checkers'}
         

        % dataTYPE name. Checkers, Words, FalseFont. 
        prf.rmName = 'Checkers'; 

        % if we want to run the model only on an roi as opposed to all gray voxels, specify path here
        prf.roiName = []; 

        % name of params file
        % p.paramsFile    = 'Stimuli/params_knkfull_multibar_blank.mat';  % Words and FalseFont
        % prf.p.paramsFile  = 'Stimuli/20150113T200620.mat';         % Checkers
        prf.p.paramsFile  = cr.bk.list_paramsFile_Checkers{subind};
        
        % image file
        % p.imFile        = 'Stimuli/images_knk_fliplr.mat';  % Words and FalseFont
        % prf.p.imFile = 'Stimuli/images_8barswithblank_fliplr.mat'; 
        prf.p.imFile      = cr.bk.list_imFile_Checkers{subind};
        
        
        % IT seems that scanNum and barScan are the same
        % scan number that with appropriate stimuli. 
        % Need to know for clipping purposes
        
        % CHECK: IS IT list_scanNum_Checkers_sessionRet
        prf.barScans      = cr.bk.list_scanNum_Checkers_sessionRet(subind);
        prf.p.scanNum     = cr.bk.list_scanNum_Checkers_sessionRet(subind);
               
    case {'Words','FalseFont'}
        % dataTYPE name. Checkers, Words, FalseFont. 
        prf.rmName = rmName; 

        % if we want to run the model only on an roi as opposed to all gray voxels, specify path here
        prf.roiName = []; 

        % name of params file
        % prf.p.paramsFile    = 'Stimuli/params_knkfull_multibar_blank.mat';  % Words and FalseFont
        % p.paramsFile  = 'Stimuli/20150113T200620.mat';         % Checkers
        prf.p.paramsFile  = cr.bk.list_paramsFile_KnK{subind};

        % image file
        % prf.p.imFile        = 'Stimuli/images_knk_fliplr.mat';  % Words and FalseFont
        % p.imFile        = 'Stimuli/images_8barswithblank_fliplr.mat'; 
        prf.p.imFile       = cr.bk.list_imFile_KnK{subind};
                
        % SEE ABOVE
        % scan number that with appropriate stimuli. 
        % Need to know for clipping purposes
        % prf.barScans  = [2]; 
        % prf.p.scanNum                  = 5;
        % CHECK: IS IT  list_scanNum_Knk_sessionRet
        prf.barScans  = cr.bk.list_scanNum_Knk_sessionRet(subind);
        prf.p.scanNum = cr.bk.list_scanNum_Knk_sessionRet(subind);
    case {'Words_English','Words_Hebrew'}
        % dataTYPE name. Checkers, Words, FalseFont. 
        prf.rmName = rmName; 

        % if we want to run the model only on an roi as opposed to all gray voxels, specify path here
        prf.roiName = []; 

        % name of params file
        % prf.p.paramsFile    = 'Stimuli/params_knkfull_multibar_blank.mat';  % Words and FalseFont
        % p.paramsFile  = 'Stimuli/20150113T200620.mat';         % Checkers
        prf.p.paramsFile  = cr.bk.list_paramsFile_KnK{subind};

        % image file
        % prf.p.imFile        = 'Stimuli/images_knk_fliplr.mat';  % Words and FalseFont
        % p.imFile        = 'Stimuli/images_8barswithblank_fliplr.mat'; 
        prf.p.imFile       = cr.bk.list_imFile_KnK{subind};
                
        % SEE ABOVE
        % scan number that with appropriate stimuli. 
        % Need to know for clipping purposes
        % prf.barScans  = [2]; 
        % prf.p.scanNum                  = 5;
        % CHECK: IS IT  list_scanNum_Knk_sessionRet
        prf.barScans  = cr.bk.list_scanNum_Knk_sessionRet(subind);
        prf.p.scanNum = cr.bk.list_scanNum_Knk_sessionRet(subind);
        
    otherwise
        error('[prfRun_params_rmName] %s not recognized',rmName)
end

