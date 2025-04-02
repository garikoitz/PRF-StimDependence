% name of params file
p.paramsFile_Knk            = 'Stimuli/params_knkfull_multibar_blank.mat';  % Words and FalseFont
p.paramsFile_Checkers       = 'Stimuli/params_checkers.mat';                % Checkers

% image file
p.imFile_Knk                = 'Stimuli/images_knk_fliplr.mat';              % Words and FalseFont
p.imFile_Checkers           = 'Stimuli/images_8barswithblank_fliplr.mat';   % Checkers

% params common to all dts
params.fliprotate       = [0 0 0]; 
params.stimType         = 'StimFromScan';
params.stimWidth        = 90;               
params.stimStart        = 0;                
params.stimDir          = 0;                
params.nCycles          = 1;               
params.nStimOnOff       = 0;                
params.nUniqueRep       = 1;                
params.nDCT             = 1;     
params.hrfType          = 'two gammas (SPM style)';
params.hrfParams        = {[1.6800 3 2.0500] [5.4000 5.2000 10.8000 7.3500 0.3500]}; 
params.imfilter         = 'binary';
params.jitterFile       = 'Stimuli/none';