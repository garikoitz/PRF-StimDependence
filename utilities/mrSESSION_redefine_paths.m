%% when switching to iMac, the global vANATOMY path references a soft link
% and iMac has problems with soft links
% So redefine and give absolute path of everyone's anatomy
% And re save the session
% use the saveSession command so that a backup will be saved
% clear all; close all; clc; 
% bookKeeping;
% cr.bk = bookKeeping(cr);
function mrSESSION_redefine_paths(cr,bk)


%% vANATOMYPATHs that have been changed
% list_sessionRet. except for subs 1 and 3
% list_sessionHebrewRet. no problems


%% modify here
% the sessions to do this for
% list_sessionRet
% list_sessionPath
% list_sessionLocPath
% list_sessionDtiQmri
% list_sessionDiffusionRun1
% list_sessionDiffusionRun2
% list_sessionHebrewRet
% list_sessionSizeRet
% list_sessionTiledLoc
% list_sessionTestRetest
% list_sessionAfq
list_session = bk.list_sessionRet; 

% subjects to do this for
% list_subInds = [1:20]% 4:38];
list_subInds = [1:length(list_session)];

%% do things
% numSessions = length(list_session); 

for ii = list_subInds
   
    %% anatomy and vista directoryes
    dirVista = list_session{ii};
    dirAnatomy = bk.list_anatomy{ii};
    
    clear GLOBALS;
    mrvCleanWorkspace;
    
    %% do things if the vista session and anatomy direcotry exist
    if exist(dirVista, 'dir') & exist(dirAnatomy, 'dir')
        
        clear dataTYPES mrSESSION vANATOMYPATH

        % move and load
        chdir(dirVista);
        mrGlobals; 
        load mrSESSION; 
        
        % the current vANATOMYPATH
        display(['old vANATOMYPATH: ' vANATOMYPATH]);
        
        % change it over and save the session
        vANATOMYPATH = fullfile(dirAnatomy, 't1.nii.gz');
        display(['new vANATOMYPATH: ' vANATOMYPATH]);        
        
        % Now do the same with the inplanes and the rest of directories
        display('Checking functionals')
        checkStructs = {'inplanes','functionals'};
        for cs = checkStructs
            if strcmp(cs{:},'inplanes')
                src      = mrSESSION.(cs{:}).inplanePath;
                srcsplit = split(src,filesep);
                oldBase  = strjoin(srcsplit(1:5,1),filesep);
                newBase  = cr.dirs.DATA;
                relLoc   = strjoin(srcsplit(end-3:end,1),filesep);
                dst      = fullfile(newBase,relLoc);
                if ~strcmp(src,dst)
                    mrSESSION.(cs{:}).inplanePath = dst;
                    display(['Changed ' cs{:}])
                end
            else
                for nm=1:length(mrSESSION.(cs{:}))
                    src      = mrSESSION.(cs{:})(nm).PfileName;
                    srcsplit = split(src,filesep);
                    oldBase  = strjoin(srcsplit(1:5,1),filesep);
                    newBase  = cr.dirs.DATA;
                    relLoc   = strjoin(srcsplit(end-3:end,1),filesep);
                    dst      = fullfile(newBase,relLoc);
                    display(['Changing ' cs{:}])
                    if ~strcmp(src,dst)
                        mrSESSION.(cs{:})(nm).PfileName = dst;
                        display(['Changed to ' dst])
                    end
                end
            end
        end
        saveSession;
        
        mrvCleanWorkspace;
        
    else
        display(['dirVista or dirAnatomy does not exist for sub ' num2str(ii)])
    end
    
end


end
