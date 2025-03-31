function Run_PRFs(cr, list_subInds)
%Run PRFs 
%   Takes the time series and runs it all again to get results in mat form
% subjects we want to do this for
% list_subInds        = [31:36 38:44];  % Hebrew
% list_subInds      = [1:20];  % Original 20
% mw (13) for Words failed, continue with the next ones for now
% list_subInds      = [18:20];
%17 and 13 failed at beginning
% list_subInds     = [1,3,4,13:20];
% list_dtNames     = {'WordSmall','WordLarge'};
% list_dtNames     = {'Checkers'};

    for subind = list_subInds

        mrvCleanWorkspace;
        % subind  = list_subInds(ns);
        subname = cr.bk.list_sub{subind};
        [~,anatName]=fileparts(cr.bk.list_anatomy{subind});
        fprintf('\nSubDetails:\nInd:%i, StrInd:%s, subname:%s, Name:%s, anatName:%s\n',...
            subind,cr.bk.list_subNumberString{subind},subname,...
            cr.bk.list_names{subind},anatName)

        % Change dir, we need to run analysis where mrSession is
        % FOR ALL
        chdir(cr.bk.list_sessionRet{subind})
        prf.dirVistacc = cr.bk.list_sessionRet{subind};
        % FOR WORD LARGE SMALL
        % chdir(cr.bk.list_sessionSizeRet{subind})
        % prf.dirVistacc = cr.bk.list_sessionSizeRet{subind};

        %% PRF analysis
        % Read the generic params for all subjects
        run(fullfile(cr.dirs.DEF,'prfrun_defaults.m'));
        cr.defaults.prfrun.params = params;
        cr.defaults.prfrun.p      = p;
        clear('params'); clear('p');
        % Read prfRun_params specific to this subject
        % run(cr.bk.list_prfParams{subind}); NOT NECESSARY

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