function sd_data = read_download_data(sd_data, cr, data_type)
%read_download_data Read or create data structure in the server
% Check if the file exists in the local dir, otherwise download from OSF
    
    for ns = 1:length(fields(sd_data))
        site_names = fields(sd_data);
        site = site_names{ns};
        
        rmroi = false; time_s = false;
        if strcmp(data_type, 'rmroiCell')
            rmroi = true; 
            fname = sd_data.(site).rmroiFname;
            fpath = fullfile(cr.dirs.local, fname);
        elseif strcmp(data_type, 'TS')
            time_s = true;
            fname = sd_data.(site).TSFname;
            fpath = fullfile(cr.dirs.local, fname);
        else
            error([data_type 'not known'])
        end
        
        % Check if exists, return it, otherwise create it
        if isfile(fpath)
            disp([fname ' exists locally, loading...'])
            if rmroi
                load(fpath,'rmroiCell')
                % Add the information to the main struct
                sd_data.(site).rmroiCell = rmroiCell;
            end
            if time_s
                load(fpath,'TS')
                % Add the information to the main struct
                sd_data.(site).TS = TS;
            end
            break
        end

        % If it is not locally available, check in OSF
        if check_download_from_osf(fpath)
            disp('rmroiCell not in local, downloading from OSF')
            if rmroi
                load(fpath,'rmroiCell')
                % Add the information to the main struct
                sd_data.(site).rmroiCell = rmroiCell;
            end
            if time_s
                load(fpath,'TS')
                % Add the information to the main struct
                sd_data.(site).TS = TS;
            end
            break
        end

        % If non of the previous options were true, create it and save it
        if rmroi
            disp('rmroiCell does not exist locally or in OSF: creating ...')
            rmroiCell = ff_rmroiCell(cr,...
                sd_data.(site).list_subInds,...
                sd_data.(site).list_roiNames,...
                sd_data.(site).list_dtNames, ...
                sd_data.(site).list_rmNames,...
                'list_path',cr.bk.list_sessionRet,...
                'latest_fFit',false, ...
                'checkYear','2022');
            % Save rmroicell, if funtcions break below, at least file is saved
            save(fpath,'rmroiCell')
            
            % Upload to OSF for the next time. 
            % upload_file_to_osf(sd_data.(site).rmroiFname, rmroiCell)
    
            % Add the information to the main struct
            sd_data.(site).rmroiCell = rmroiCell;
        end
        if time_s
            disp('TS does not exist locally or in OSF: creating ...')
            TS = ff_rmroiCell(cr,...
                        sd_data.(site).list_subInds,...
                        sd_data.(site).list_roiNames,...
                        sd_data.(site).list_dtNames, ...
                        sd_data.(site).list_rmNames,...
                        'list_path',cr.bk.list_sessionRet,...
                        'latest_fFit',false, ...
                        'checkYear','2022', ...
                        'data_type', data_type);
            % Save rmroicell, if funtcions break below, at least file is saved
            save(fpath,'TS')
            
            % Upload to OSF for the next time. 
            % upload_file_to_osf(sd_data.(site).rmroiFname, rmroiCell)
    
            % Add the information to the main struct
            sd_data.(site).TS = TS;
        end
    end
end
