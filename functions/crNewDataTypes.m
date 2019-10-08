function cr = crNewDataTypes(cr, opt, varargin)
% function that will create a new dataTYPES that is the average of existing scans
% This function is based on the script s_tSeriesAverageAndXform_individualRuns.m
% by Rosemary, created in coverageReading

% close all; clear all; clc; 
% bookKeeping; 
% dirVista = '/sni-storage/wandell/data/reading_prf/heb_pilot09/RetAndHebrewLoc';

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


% TODO: clean this duplicates of filenames

% Add other lists coming from bookKeeping
% The following parameters have been set in bk = bookKeeping():
% session list. see bookKeeping
opt.list_path            = cr.bk.list_sessionRet; 
% list of checker scan number, corresponding to session list
opt.list_numCheckers     = cr.bk.list_scanNum_Checkers_sessionRet; 
% list of knk scan number, corresponding to session list
opt.list_numKnk          = cr.bk.list_scanNum_Knk_sessionRet;  

%% Don't loop over all subjects. Do that from the main script
% Every subject specific value has to be written in bookKeeping and should be
% extracted from there. 

% Specify session path. usually this script is saved right there 
cr.dirVista = opt.list_path{opt.list_subInds};
chdir(cr.dirVista);

% create params structure ...
params = mrInitDefaultParams; 

%% 

vwI = initHiddenInplane; 
vwG = initHiddenGray; 

%%
% check that dtAssignments matches the number of scans we have
% set to a motion corrected dt, which should have total number of sans
vwI         = viewSet(vwI,'curdt',opt.dtToAverage);
numScans    = viewGet(vwI,'numScans');
if numScans ~= length(opt.dtAssignments)
    error('Scan length mismatch!')
end

for ii = 1:length(opt.dtsToCreate)
    
    % the datatype that we will create the new dt from.
    vwI  = viewSet(vwI,'curdt',opt.dtToAverage);
    
    % datatype name
    newDtName = opt.dtsToCreate{ii}; 
    
    % which scans to average
    scansToAvg = find(opt.dtAssignments == ii); 
    
    % average
    vwI  =  averageTSeries(vwI, scansToAvg, newDtName);
    
    % set inplane and gray newDtName
    vwI  = viewSet(vwI,'curdt',newDtName);
    vwG  = viewSet(vwG,'curdt', newDtName); 
    
    % do the xform
    ip2volTSeries(vwI,vwG,0,'linear'); 
    
end

end
