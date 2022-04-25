% This scripts will download the dicoms from Flywheel, and then convert them to
% niftis using heudiconv and then run fmriprep, so that we can reanalyze
% everything using prfprepare and prfanalyze-vista

% Variables to edit
serverName     = 'stanfordlabs';
collectionName = 'VWFA_FOV';
showSessions   = false;


%% 1.- CONNECT AND READ SESSIONS (server & collection)

st = scitran(serverName);
st.verify

% Connect to the collection, verify it and show the number of sessions for verification
% FC: obtain collection ID from the collection name
collectionID = '';
collections  = st.fw.getAllCollections();
for nc=1:length(collections)
    if strcmp(collections{nc}.label, collectionName)
        collectionID = collections{nc}.id;
    end
end

if isempty(collectionID)
    error(sprintf('Collection %s could not be found on the server %s (verify permissions or the collection name).', collectionName, serverName))
else
    thisCollection        = st.fw.getCollection(collectionID);
    sessionsInCollection  = st.fw.getCollectionSessions(idGet(thisCollection));
    fprintf('There are %i sessions in the collection %s (server %s).\n', length(sessionsInCollection), collectionName, serverName)
    if showSessions
        for ns=1:length(sessionsInCollection)
            thisSession = st.fw.getSession(idGet(sessionsInCollection{ns}));
            % Get info for the project the session belong to
            thisProject = st.fw.getProject(thisSession.project);
            fprintf('(%d) %s >> %s (%s)\n', ns, thisProject.label, thisSession.subject.code, thisSession.label)
        end
    end
end

% Read the bookkeeping file from Rosemary too
cr         = struct();
cr.codeDir = crRP;

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
% addpath(genpath(fullfile(cr.dirs.ANALYSIS,'matlabfiles')));

% Rosemary relied on this file that contains most of the subjects and other
% lists. Make it work with relative paths and store it in each project repository
% This file was used as well for: 
% - copying files to a new location
% - editing mrSession to reflect the file changes
cr.bk = bookKeeping(cr);

%% 2.- Create table with files within each session
% Find the subjects we are interested in FW
list_subInds        = [1:20];  % Original 20
% list_subInds        = [31:36 38:44];  % Hebrew

names20 = {};
ii = 0;
for subind = list_subInds
    ii = ii + 1;
    subname = cr.bk.list_sub{subind};
    [~,anatName]=fileparts(cr.bk.list_anatomy{subind});
    fprintf('\nSubDetails:\nInd:%i, StrInd:%s, subname:%s, Name:%s, anatName:%s\n',...
        subind,cr.bk.list_subNumberString{subind},subname,...
        cr.bk.list_names{subind},anatName)
    names20{ii} = lower(cr.bk.list_names{subind});
end





namesurs = {};
for ns=1:length(sessionsInCollection)
    thisSession = st.fw.getSession(idGet(sessionsInCollection{ns}));
    % Get info for the project the session belong to
    thisProject  = st.fw.getProject(thisSession.project);
    projectLabel = thisProject.label;
    subCode      = thisSession.subject.code;
    subj         = st.fw.getSubject(idGet(thisSession.subject));
%     sesCode      = thisSession.label;
%     acqus        = st.list('acquisition', idGet(thisSession));
%     for na = 1:length(acqus)
%         thisAcqu = st.fw.getAcquisition(idGet(acqus{na}));
%         thisAcqu.label
%     end
    namesur = [lower(subj.firstname) ' ' lower(subj.lastname)];
    if ismember(namesur, names20)  
        fprintf('\nSubDetails: Name:%s \n',namesur)
    end
    namesurs{ns} = namesur;
end

length(unique(namesurs))





for ns = 1: length(sessionsInCollection)
sessionsInCollection{1}.subject.info
a = thisSession.acquisitions()
a{1}.files{1}
