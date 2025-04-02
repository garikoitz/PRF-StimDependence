function is_file_osf = check_download_from_osf(rmroi_fpath)
%check_download_from_osf Check if file is in osf, if it is, download it. 
    is_file_osf = false;

    [fpath, fname, fext] = fileparts(rmroi_fpath);
    switch fname
        case 'rmroicell_subInds-1to20_dtNames-cb-w-ff_fits-new_LeftRightROIs_2023'
            url = 'https://osf.io/download/y7wp6/';
        case 'rmroicell_subInds-31to36-38to44_dtNames-ALL-LeftRight_fits-new_2023'
            url = 'https://osf.io//download/b35qe/';
        case 'PRFstim_sd_data_31march2025'
            url = 'https://osf.io/download/9dawj';
        case 'TS_CNI_subInds-1-20_dtNames-CB-WE-FF_ROI-all'
            url = 'https://osf.io/download/vzbuj';
        case 'TS_ISRAEL_subInds-31-36-38-44_dtNames-CB-WE-WH_ROI-all'
            url = 'https://osf.io/download/aq62c';
        case 'rmroicell_CNI_subInds-1-20_dtNames-CB-WE-FF_ROI-all_fits-fFit_2025'
            url = 'https://osf.io/download/9yz2b';
        case 'rmroicell_ISRAEL_subInds-31-36-38-44_dtNames-CB-WE-WH_ROI-all_fits-fFit_2025'
            url = 'https://osf.io/download/a2vwz';
        otherwise
            error([fname fext ' file not known'])
    end

    websave(rmroi_fpath, url);
    is_file_osf = true;
 
end



