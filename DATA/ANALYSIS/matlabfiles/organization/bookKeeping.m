function bk = bookKeeping(cr)

% curdir = pwd;
dataDir     = cr.dirs.DATA;
defaultsDir = cr.dirs.DEF;
checkCopy   = false;
redefinePath= false;
export2fw   = false;  % do it with michael perry to upload as analysis some of it



%% Move the old data to the new location where we control de script
% We will delete this at the end, use this to copy files. 
% Maintain the code somewhere for reference.
oldbaseDir = '/biac4/wandell/data';
olddataDir = fullfile(oldbaseDir, 'reading_prf');
oldanatomy    = '/biac4/wandell/data/reading_prf/anatomy';

%% bookkeeping for retinotopy and words analysis
% chdir('/share/wandell/data/reading_prf/coverageReading/organization')
% chdir('/share/wandell/data/reading_prf/coverageReading/organization')
% chdir(fullfile(crRootPath, 'organization'));

%% name of sessions, abbreviations, and their paths
% subject initials
bk.list_sub = {
    'jg'        % 1
    'ad'        % 2
    'cc'        % 3
    'jw'        % 4
    'rs'        % 5
    'sg'        % 6
    'th'        % 7
    'pv'        % 8
    'sl'        % 9
    'jv'        % 10
    'dl'        % 11
    'ak'        % 12
    'mw'        % 13
    'gt'        % 14
    'ws'        % 15
    'ol'        % 16
    'tl'        % 17
    'mv'        % 18
    'vm'        % 19
    'ab'        % 20
    'bw'        % 21
    'dys_ab'    % 22
    'heb_ag'    % 23
    'heb_aa'    % 24
    'heb_ls'    % 25
    'heb_toba'  % 26
    'heb_maya'  % 27
    'heb_yama'  % 28
    'heb_blta'  % 29
    'heb_mibe'  % 30
    'heb_avbe'  % 31
    'heb_nitr'  % 32
    'heb_gilo'  % 33
    'heb_eisa'  % 34
    'heb_dael'  % 35
    'heb_taay'  % 36
    'heb_maaf'  % 37
    'heb_brne'  % 38
    'heb_nagr'  % 39
    'heb_avhi'  % 40
    'heb_avar'  % 41
    'heb_nihe'  % 42
    'heb_noco'  % 43
    'heb_shsh'  % 44
    };


bk.list_subNumberString = {
    '01'
    '02'
    '03'
    '04'
    '05'
    '06'
    '07'
    '08'
    '09'
    '10'
    '11'
    '12'
    '13'
    '14'
    '15'
    '16'
    '17'
    '18'
    '19'
    '20'
    '21'
    '22'
    '23'
    '24'
    '25'
    '26'
    '27'
    '28'
    '29'
    '30'
    '31'
    '32'
    '33'
    '34'
    '35'
    '36'
    '37'
    '38'
    '39'
    '40'
    '41'
    '42'
    '43'
    '44'
    };

bk.list_names = {
    'Jesse Gomez'        % 1
    'Aaron Dames'        % 2
    'Cat Camacho'        % 3
    'Joe Wexler'         % 4
    'Rose Schneider'     % 5
    'Sarah Gleberman'    % 6
    'Tyler Hughes'       % 7
    'Priyanka Vij'       % 8
    'Stefanie Lim'       % 9
    'Janet Veil'         % 10
    'Destiny Lopez'      % 11
    'Anna Khazenzon'     % 12
    'Michael Waskom'     % 13
    'Gerald Tiu'         % 14
    'Wendy Sato'         % 15
    'Oliver Leung'       % 16
    'Trisha Lian'        % 17
    'Michael Vitelli'    % 18
    'Vanessa Martinez'   % 19
    'Annette Bugno'      % 20
    'Brian Wandell'      % 21
    'Andy Brainard'      % 22   Dyslexic
    'Amit Goodman'       % 23   Hebrew
    'Avital Ayzenshtat'  % 24   Hebrew
    'Linoi Shambiro'     % 25   Hebrew
    'To Ba'              % 26   Hebrew
    'Ma Ya'              % 27   Hebrew
    'Ya Ma'              % 28   Hebrew
    'Bl Ta'              % 29   Hebrew
    'Michal Ben-Shachar' % 30   Hebrew
    'Av Be'              % 31   Hebrew
    'Ni Tr'              % 32   Hebrew
    'Gi Lo'              % 33   Hebrew
    'Ei Sa'              % 34   Hebrew
    'Da El'              % 35   Hebrew
    'Ta Ay'              % 36   Hebrew
    'Ma Af'              % 37   Hebrew
    'Br Ne'              % 38   Hebrew
    'Na Gr'              % 39   Hebrew
    'Av Hi'              % 40   Hebrew
    'Av Ar'              % 41   Hebrew
    'Ni He'              % 42   Hebrew
    'No Co'              % 43   Hebrew
    'Sh Sh'              % 44   Hebrew
    };

% directory with ret. Checkers. English Words
bk.list_sessionRet = {
    fullfile(dataDir,'jg','20150113_1947');       % jg
    fullfile(dataDir,'ad','20150120_ret');        % ad
    fullfile(dataDir,'cc','20150122_ret');        % cc
    fullfile(dataDir,'jw','20150124_ret');        % jw
    fullfile(dataDir,'rs','20150201_ret');        % rs
    fullfile(dataDir,'sg','20150131_ret');        % sg
    fullfile(dataDir,'th','20150205_ret');        % th
    fullfile(dataDir,'pv','20150503_ret');        % pv
    fullfile(dataDir,'sl','20150507_ret');        % sl
    fullfile(dataDir,'jv','20150509_ret');        % jv
    fullfile(dataDir,'dl','20150519_ret');        % dl
    fullfile(dataDir,'ak','20150106_1802');       % ak
    fullfile(dataDir,'mw','tiledLoc_sizeRet');    % mw
    fullfile(dataDir,'gt','tiledLoc_sizeRet');    % gt
    fullfile(dataDir,'ws','tiledLoc_sizeRet');    % ws
    fullfile(dataDir,'ol','tiledLoc_sizeRet');    % ol
    fullfile(dataDir,'tl','Localizer_sizeRet');   % tl% the stim file is wrong, check it
    fullfile(dataDir,'mv','tiledLoc_sizeRet');    % mv
    fullfile(dataDir,'vm','tiledLoc_sizeRet');    % vm
    fullfile(dataDir,'ab','tiledLoc_sizeRet');    % ab
    fullfile(dataDir,'bw','tiledLoc_sizeRet');    % bw
    fullfile(dataDir,'dys_ab','20150430_ret');    % dys_ab
    fullfile(dataDir,'heb_pilot01','Analyze_pseudoInplane');      % heb_ag
    fullfile(dataDir,'heb_pilot02','RetAndLoc');                  % heb_aa
    fullfile(dataDir,'heb_pilot03','RetAndLoc');                  % heb_ls
    fullfile(dataDir,'heb_pilot04','RetAndLoc_noXform');          % heb_toba
    fullfile(dataDir,'heb_pilot05','RetAndLoc');                  % heb_maya
    fullfile(dataDir,'heb_pilot06','Ret_OTSPrescription_152Vol'); % heb_yama
    fullfile(dataDir,'heb_pilot07','RetAndHebrewLoc');            % heb_blta
    fullfile(dataDir,'heb_pilot08','RetAndHebrewLoc');            % heb_mibe
    fullfile(dataDir,'heb_pilot09','RetAndHebrewLoc');             % 31
    fullfile(dataDir,'heb_pilot10','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot11','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot12','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot13','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot14','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot15','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot16','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot17','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot18','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot19','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot20','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot21','RetAndHebrewLoc');
    fullfile(dataDir,'heb_pilot22','RetAndHebrewLoc');
};

if checkCopy,checkAndCopy(bk.list_sessionRet, bk.list_sessionRet),end


% Added by gari
bk.list_prfParams = {
    fullfile(dataDir,'jg','20150113_1947','params','prfRun_params.m');       % jg
    fullfile(dataDir,'ad','20150120_ret','params','prfRun_params.m');        % ad
    fullfile(dataDir,'cc','20150122_ret','params','prfRun_params.m');        % cc
    fullfile(dataDir,'jw','20150124_ret','params','prfRun_params.m');        % jw
    fullfile(dataDir,'rs','20150201_ret','params','prfRun_params.m');        % rs
    fullfile(dataDir,'sg','20150131_ret','params','prfRun_params.m');        % sg
    fullfile(dataDir,'th','20150205_ret','params','prfRun_params.m');        % th
    fullfile(dataDir,'pv','20150503_ret','params','prfRun_params.m');        % pv
    fullfile(dataDir,'sl','20150507_ret','params','prfRun_params.m');        % sl
    fullfile(dataDir,'jv','20150509_ret','params','prfRun_params.m');        % jv
    fullfile(dataDir,'dl','20150519_ret','params','prfRun_params.m');        % dl
    fullfile(dataDir,'ak','20150106_1802','params','prfRun_params.m');       % ak
    fullfile(dataDir,'mw','tiledLoc_sizeRet','params','prfRun_params.m');    % mw
    fullfile(dataDir,'gt','tiledLoc_sizeRet','params','prfRun_params.m');    % gt
    fullfile(dataDir,'ws','tiledLoc_sizeRet','params','prfRun_params.m');    % ws
    fullfile(dataDir,'ol','tiledLoc_sizeRet','params','prfRun_params.m');    % ol
    fullfile(dataDir,'tl','Localizer_sizeRet','params','prfRun_params.m');   % tl
    fullfile(dataDir,'mv','tiledLoc_sizeRet','params','prfRun_params.m');    % mv
    fullfile(dataDir,'vm','tiledLoc_sizeRet','params','prfRun_params.m');    % vm
    fullfile(dataDir,'ab','tiledLoc_sizeRet','params','prfRun_params.m');    % ab
    fullfile(dataDir,'bw','tiledLoc_sizeRet','params','prfRun_params.m');    % bw
    fullfile(dataDir,'dys_ab','20150430_ret','params','prfRun_params.m');    % dys_ab
    fullfile(dataDir,'heb_pilot01','Analyze_pseudoInplane','params','prfRun_params.m');      % heb_ag
    fullfile(dataDir,'heb_pilot02','RetAndLoc','params','prfRun_params.m');                  % heb_aa
    fullfile(dataDir,'heb_pilot03','RetAndLoc','params','prfRun_params.m');                  % heb_ls
    fullfile(dataDir,'heb_pilot04','RetAndLoc_noXform','params','prfRun_params.m');          % heb_toba
    fullfile(dataDir,'heb_pilot05','RetAndLoc','params','prfRun_params.m');                  % heb_maya
    fullfile(dataDir,'heb_pilot06','Ret_OTSPrescription_152Vol','params','prfRun_params.m'); % heb_yama
    fullfile(dataDir,'heb_pilot07','RetAndHebrewLoc','params','prfRun_params.m');            % heb_blta
    fullfile(dataDir,'heb_pilot08','RetAndHebrewLoc','params','prfRun_params.m');            % heb_mibe
    fullfile(dataDir,'heb_pilot09','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot10','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot11','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot12','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot13','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot14','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot15','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot16','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot17','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot18','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot19','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot20','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot21','RetAndHebrewLoc','params','prfRun_params.m');
    fullfile(dataDir,'heb_pilot22','RetAndHebrewLoc','params','prfRun_params.m');
};

% Added by gari
bk.list_coverageFigure_defaults = {
    fullfile(dataDir,'jg','20150113_1947','params','coverageFigure_defaults.m');       % jg
    fullfile(dataDir,'ad','20150120_ret','params','coverageFigure_defaults.m');        % ad
    fullfile(dataDir,'cc','20150122_ret','params','coverageFigure_defaults.m');        % cc
    fullfile(dataDir,'jw','20150124_ret','params','coverageFigure_defaults.m');        % jw
    fullfile(dataDir,'rs','20150201_ret','params','coverageFigure_defaults.m');        % rs
    fullfile(dataDir,'sg','20150131_ret','params','coverageFigure_defaults.m');        % sg
    fullfile(dataDir,'th','20150205_ret','params','coverageFigure_defaults.m');        % th
    fullfile(dataDir,'pv','20150503_ret','params','coverageFigure_defaults.m');        % pv
    fullfile(dataDir,'sl','20150507_ret','params','coverageFigure_defaults.m');        % sl
    fullfile(dataDir,'jv','20150509_ret','params','coverageFigure_defaults.m');        % jv
    fullfile(dataDir,'dl','20150519_ret','params','coverageFigure_defaults.m');        % dl
    fullfile(dataDir,'ak','20150106_1802','params','coverageFigure_defaults.m');       % ak
    fullfile(dataDir,'mw','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % mw
    fullfile(dataDir,'gt','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % gt
    fullfile(dataDir,'ws','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % ws
    fullfile(dataDir,'ol','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % ol
    fullfile(dataDir,'tl','Localizer_sizeRet','params','coverageFigure_defaults.m');   % tl
    fullfile(dataDir,'mv','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % mv
    fullfile(dataDir,'vm','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % vm
    fullfile(dataDir,'ab','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % ab
    fullfile(dataDir,'bw','tiledLoc_sizeRet','params','coverageFigure_defaults.m');    % bw
    fullfile(dataDir,'dys_ab','20150430_ret','params','coverageFigure_defaults.m');    % dys_ab
    fullfile(dataDir,'heb_pilot01','Analyze_pseudoInplane','params','coverageFigure_defaults.m');      % heb_ag
    fullfile(dataDir,'heb_pilot02','RetAndLoc','params','coverageFigure_defaults.m');                  % heb_aa
    fullfile(dataDir,'heb_pilot03','RetAndLoc','params','coverageFigure_defaults.m');                  % heb_ls
    fullfile(dataDir,'heb_pilot04','RetAndLoc_noXform','params','coverageFigure_defaults.m');          % heb_toba
    fullfile(dataDir,'heb_pilot05','RetAndLoc','params','coverageFigure_defaults.m');                  % heb_maya
    fullfile(dataDir,'heb_pilot06','Ret_OTSPrescription_152Vol','params','coverageFigure_defaults.m'); % heb_yama
    fullfile(dataDir,'heb_pilot07','RetAndHebrewLoc','params','coverageFigure_defaults.m');            % heb_blta
    fullfile(dataDir,'heb_pilot08','RetAndHebrewLoc','params','coverageFigure_defaults.m');            % heb_mibe
    fullfile(dataDir,'heb_pilot09','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot10','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot11','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot12','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot13','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot14','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot15','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot16','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot17','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot18','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot19','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot20','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot21','RetAndHebrewLoc','params','coverageFigure_defaults.m');
    fullfile(dataDir,'heb_pilot22','RetAndHebrewLoc','params','coverageFigure_defaults.m');
};

% directory with checkers, words, falsefont retinotopy
bk.list_sessionPath = {
    % '/share/wandell/data/reading_prf/rosemary/20141026_1148';     % rl
    fullfile(dataDir,'jg','20150113_1947');       % jg
    fullfile(dataDir,'ad','20150120_ret');        % ad
    fullfile(dataDir,'cc','20150122_ret');        % cc
    fullfile(dataDir,'jw','20150124_ret');        % jw
    fullfile(dataDir,'rs','20150201_ret');        % rs
    fullfile(dataDir,'sg','20150131_ret');        % sg
    fullfile(dataDir,'th','20150205_ret');        % th
    fullfile(dataDir,'pv','20150503_ret');        % pv
    fullfile(dataDir,'sl','20150507_ret');        % sl
    fullfile(dataDir,'jv','20150509_ret');        % jv
    fullfile(dataDir,'dl','20150519_ret');        % dl
    fullfile(dataDir,'ak','20150106_1802')        % ak
    'do not have'                                                   % mw
    'do not have'                                                   % gt
    'do not have'                                                   % ws
    'do not have'                                                   % ol
    'do not have'                                                   % tl
    'do not have'                                                   % mv
    'do not have'                                                   % vm
    'do not have'                                                   % ab
    'do not have'                                                   % bw
    fullfile(dataDir,'dys_ab','20150430_ret');    % dys_ab
    'do not have'                                                   % heb_ag
    'do not have'                                                   % heb_aa
    'do not have'                                                   % heb_ls
    'do not have'                                                   % heb_toba
    'do not have'                                                   % heb_maya
    'do not have'                                                   % heb_yama
    'do not have'                                                   % heb_blta
    'do not have'                                                   % heb_mibe
    'do not have'                                                   % heb_avbe
    'do not have'                                                   % heb_nitr
    'do not have'                                                   % heb_gilo
    'do not have'                                                   % heb_eisa
    'do not have'                                                   % heb_dael
    'do not have'                                                   % heb_taay
    'do not have'                                                   % heb_maaf
    'do not have'                                                   % heb_brne
    'do not have'                                                   % heb_nagr
    'do not have'                                                   % heb_avhi
    'do not have'                                                   % heb_avar
    'do not have'                                                   % heb_nihe
    'do not have'                                                   % heb_noco
    'do not have'                                                   % heb_shsh
    };
if checkCopy,checkAndCopy(bk.list_sessionPath, bk.list_sessionRet),end


% directory with the smallFOV localizer mrSESSION
bk.list_sessionLocPath = {
    % '/Volumes/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Localizer/data/RL22_05312014'; % rk
    % '/Volumes/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Localizer/data/JG24_01062014';   % jg, move it
    fullfile(dataDir,'jg','JG24_01062014');                              % jg
    fullfile(dataDir,'ad','20150120_localizer');                              % ad
    fullfile(dataDir,'cc','20150122_loc');                                    % cc
    fullfile(dataDir,'jw','20150126_loc');                                    % jw
    fullfile(dataDir,'rs','20150201_loc');                                    % rs
    fullfile(dataDir,'sg','20150131_loc');                                    % sg
    fullfile(dataDir,'th','20150205_loc');                                    % th
    fullfile(dataDir,'pv','20150503_loc');                                    % pv
    fullfile(dataDir,'sl','20150507_loc');                                    % sl
    fullfile(dataDir,'jv','20150509_loc');                                    % jv
    fullfile(dataDir,'dl','20150519_loc');                                    % dl
    fullfile(dataDir,'ak','20150102_1039')                                    % ak
    'do not have';                                                                              % mw
    'do not have';                                                                              % gt
    'do not have';                                                                              % ws
    'do not have';                                                                              % ol
    'do not have';                                                                              % tl
    'do not have';                                                                              % mv
    'do not have';                                                                              % vm
    'do not have';                                                                              % ab
    'do not have';                                                                              % bw
    fullfile(dataDir,'dys_ab','20150430_loc');                                      % dys_ab
    'do not have';                                                                               % heb_ag
    'do not have';                                                                               % heb_aa
    'do not have';                                                                               % heb_ls
    'do not have';                                                                               % heb_toba
    'do not have';                                                                               % heb_maya
    'do not have';                                                                               % heb_yama
    'do not have';                                                                               % heb_blta
    'do not have';                                                                               % heb_mibe
    'do not have';                                                                               % heb_avbe
    'do not have';                                                                               % heb_nitr
    'do not have';                                                                               % heb_gilo
    'do not have';                                                                               % heb_eisa
    'do not have';                                                                               % heb_dael
    'do not have';                                                                               % heb_taay
    'do not have';                                                                               % heb_maaf
    'do not have';                                                                               % heb_brne
    'do not have';                                                                               % heb_nagr
    'do not have';                                                                               % heb_avhi
    'do not have';                                                                            % heb_avar
    'do not have';                                                                                % heb_nihe
    'do not have';                                                                                % heb_noco
    'do not have';                                                                                % heb_shsh
    };

if checkCopy,checkAndCopy(bk.list_sessionLocPath, bk.list_sessionRet),end

% {
% for backwawrds compatibility this list is kept. list_sessionDiffusionRun1
% and list_sessionDiffusionRun2 are better for code readbility
% directory with the dti
bk.list_sessionDtiQmri = {
    'do not have';                                                % jg
    fullfile(dataDir,'ad','diffusion');               % ad
    fullfile(dataDir,'cc','diffusion');               % cc
    fullfile(dataDir,'jw','diffusion');               % jw
    fullfile(dataDir,'rs','diffusion');               % rs
    fullfile(dataDir,'sg','diffusion');               % sg
    fullfile(dataDir,'th','diffusion');               % th
    fullfile(dataDir,'pv','diffusion');               % pv
    fullfile(dataDir,'sl','diffusion');               % sl
    fullfile(dataDir,'jv','diffusion');               % jv
    'do not have';                                                % dl
    'do not have';                                                % ak
    fullfile(dataDir,'mw','diffusion');               % mw
    fullfile(dataDir,'gt','diffusion');               % gt
    fullfile(dataDir,'ws','diffusion');               % ws
    fullfile(dataDir,'ol','diffusion');               % ol
    fullfile(dataDir,'tl','diffusion');               % tl
    fullfile(dataDir,'mv','diffusion');               % mv
    'do not have';                                                  % vm
    'do not have';                                                  % ab
    'do not have';                                                  % bw
    fullfile(dataDir,'dys_ab','diffusion');             % dys_ab
    'do not have';                                                  % heb_ag
    'do not have';                                                  % heb_aa
    'do not have';                                                  % heb_ls
    'do not have';                                                  % heb_toba
    'do not have';                                                  % heb_maya
    'do not have';                                                  % heb_yama
    'do not have';                                                  % heb_blta
    'do not have';                                                  % heb_mibe
    'do not have';                                                  % heb_avbe
    'do not have';                                                  % heb_nitr
    'do not have';                                                  % heb_gilo
    'do not have';                                                  % heb_eisa
    'do not have';                                                  % heb_dael
    'do not have';                                                  % heb_taay
    'do not have';                                                  % heb_maaf
    'do not have';                                                  % heb_brne
    'do not have';                                                  % heb_avhi
    'do not have';                                                  % heb_nagr
    'do not have'                                                   % heb_avar
    'do not have'                                                   % heb_nihe
    'do not have'                                                   % heb_noco
    'do not have'                                                   % heb_shsh
    };
if checkCopy,checkAndCopy(bk.list_sessionDtiQmri, bk.list_sessionRet),end

% directories with RUN 1 of diffusion data
bk.list_sessionDiffusionRun1 = {
    'do not have';                                                % jg
    fullfile(dataDir,'ad','diffusion');               % ad
    fullfile(dataDir,'cc','diffusion');               % cc
    fullfile(dataDir,'jw','diffusion');               % jw
    fullfile(dataDir,'rs','diffusion');               % rs
    fullfile(dataDir,'sg','diffusion');               % sg
    fullfile(dataDir,'th','diffusion');               % th
    fullfile(dataDir,'pv','diffusion');               % pv
    fullfile(dataDir,'sl','diffusion');               % sl
    fullfile(dataDir,'jv','diffusion');               % jv
    'do not have';                                                % dl
    'do not have';                                                % ak
    fullfile(dataDir,'mw','diffusion');               % mw
    fullfile(dataDir,'gt','diffusion');               % gt
    fullfile(dataDir,'ws','diffusion');               % ws
    fullfile(dataDir,'ol','diffusion');               % ol
    fullfile(dataDir,'tl','diffusion');               % tl
    fullfile(dataDir,'mv','diffusion');               % mv
    'do not have';                                                % vm
    'do not have';                                                % ab
    'do not have';                                                % bw
    fullfile(dataDir,'dys_ab','diffusion');           % dys_ab
    'do not have';                                                  % heb_ag
    'do not have';                                                  % heb_aa
    'do not have';                                                  % heb_ls
    'do not have';                                                  % heb_toba
    'do not have';                                                  % heb_maya
    'do not have';                                                  % heb_yama
    'do not have';                                                  % heb_blta
    'do not have';                                                  % heb_mibe
    'do not have';                                                  % heb_avbe
    'do not have';                                                  % heb_nitr
    'do not have';                                                  % heb_gilo
    'do not have';                                                  % heb_eisa
    'do not have';                                                  % heb_dael
    'do not have';                                                  % heb_taay
    'do not have';                                                  % heb_maaf
    'do not have';                                                  % heb_brne
    'do not have';                                                  % heb_nagr
    'do not have';                                                  % heb_avhi
    'do not have'                                                   % heb_avar
    'do not have'                                                   % heb_nihe
    'do not have'                                                   % heb_noco
    'do not have'                                                   % heb_shsh
    };
if checkCopy,checkAndCopy(bk.list_sessionDiffusionRun1,bk.list_sessionDtiQmri),end

% directories with RUN 2 of diffusion data
bk.list_sessionDiffusionRun2 = {
    'do not have';                                                % jg
    fullfile(dataDir,'ad','diffusionRun2');           % ad
    fullfile(dataDir,'cc','diffusionRun2');           % cc
    fullfile(dataDir,'jw','diffusionRun2');           % jw
    fullfile(dataDir,'rs','diffusionRun2');           % rs
    fullfile(dataDir,'sg','diffusionRun2');           % sg
    fullfile(dataDir,'th','diffusionRun2');           % th
    fullfile(dataDir,'pv','diffusionRun2');           % pv
    fullfile(dataDir,'sl','diffusionRun2');           % sl
    fullfile(dataDir,'jv','diffusionRun2');           % jv
    'do not have';                                                  % dl
    'do not have';                                                  % ak
    fullfile(dataDir,'mw','diffusionRun2');           % mw
    fullfile(dataDir,'gt','diffusionRun2');           % gt
    fullfile(dataDir,'ws','diffusionRun2');           % ws
    fullfile(dataDir,'ol','diffusionRun2');           % ol
    fullfile(dataDir,'tl','diffusionRun2');           % tl
    fullfile(dataDir,'mv','diffusionRun2');           % mv
    'do not have';                                                  % vm
    'do not have';                                                  % ab
    'do not have';                                                  % bw
    fullfile(dataDir,'dys_ab','diffusionRun2');       % dys_ab
    'do not have';                                                  % heb_ag
    'do not have';                                                  % heb_aa
    'do not have';                                                  % heb_ls
    'do not have';                                                  % heb_toba
    'do not have';                                                  % heb_maya
    'do not have';                                                  % heb_yama
    'do not have';                                                  % heb_blta
    'do not have';                                                  % heb_mibe
    'do not have';                                                  % heb_avbe
    'do not have';                                                  % heb_nitr
    'do not have';                                                  % heb_gilo
    'do not have';                                                  % heb_eisa
    'do not have';                                                  % heb_dael
    'do not have';                                                  % heb_taay
    'do not have';                                                  % heb_maaf
    'do not have';                                                  % heb_brne
    'do not have';                                                  % heb_nagr
    'do not have';                                                  % heb_avhi
    'do not have'                                                   % heb_avar
    'do not have'                                                   % heb_nihe
    'do not have'                                                   % heb_noco
    'do not have'                                                   % heb_shsh
    };
if checkCopy,checkAndCopy(bk.list_sessionDiffusionRun2, bk.list_sessionDiffusionRun2),end

% session with hebrew retinotopy. 4.5 degree radius
bk.list_sessionHebrewRet = {
    '';                  % jg
    '';                  % ad
    '';                  % cc
    '';                  % jw
    '';                  % rs
    '';                  % sg
    '';                  % th
    '';                  % pv
    '';                  % sl
    '';                  % jv
    '';                  % dl
    '';                  % ak
    '';                  % mw
    '';                  % gt
    '';                  % ws
    '';                  % ol
    '';                  % tl
    '';                  % mv
    '';                  % vm
    '';                  % ab
    '';                  % bw
    '';                  % dys_ab
    '';                  % heb_ag
    '';                  % heb_aa
    '';                  % heb_ls
    '';                  % heb_toba
    '';                  % heb_maya
    fullfile(dataDir,'heb_pilot06','Ret_OTSPrescription_152Vol'); % heb_yama
    fullfile(dataDir,'heb_pilot07','RetAndHebrewLoc');            % heb_blta
    fullfile(dataDir,'heb_pilot08','RetAndHebrewLoc');            % heb_mibe
    fullfile(dataDir,'heb_pilot09','RetAndHebrewLoc');            % heb_avbe
    fullfile(dataDir,'heb_pilot10','RetAndHebrewLoc');            % heb_nitr
    fullfile(dataDir,'heb_pilot11','RetAndHebrewLoc');            % heb_gilo
    fullfile(dataDir,'heb_pilot12','RetAndHebrewLoc');            % heb_eisa
    fullfile(dataDir,'heb_pilot13','RetAndHebrewLoc');            % heb_dael
    fullfile(dataDir,'heb_pilot14','RetAndHebrewLoc');            % heb_taay
    fullfile(dataDir,'heb_pilot15','RetAndHebrewLoc');            % heb_maaf
    fullfile(dataDir,'heb_pilot16','RetAndHebrewLoc');            % heb_brne
    fullfile(dataDir,'heb_pilot17','RetAndHebrewLoc');            % heb_nagr
    fullfile(dataDir,'heb_pilot18','RetAndHebrewLoc');            % heb_avhi
    fullfile(dataDir,'heb_pilot19','RetAndHebrewLoc');            % heb_avar
    fullfile(dataDir,'heb_pilot20','RetAndHebrewLoc');            % heb_nihe
    fullfile(dataDir,'heb_pilot21','RetAndHebrewLoc');            % heb_noco
    fullfile(dataDir,'heb_pilot22','RetAndHebrewLoc');            % heb_shsh
    };

if checkCopy,checkAndCopy(bk.list_sessionHebrewRet, bk.list_sessionRet),end
% anatomy directory. With mrVista, this information is in vANATOMYPATH, but
% when intializing for dti data it is easier if this info is just specified
% here
bk.list_anatomy = {
    fullfile(dataDir,'anatomy','gomez');        % jg
    fullfile(dataDir,'anatomy','dames');        % ad
    fullfile(dataDir,'anatomy','camacho');      % cc
    fullfile(dataDir,'anatomy','wexler');       % jw
    fullfile(dataDir,'anatomy','schneider');    % rs
    fullfile(dataDir,'anatomy','gleberman');    % sg
    fullfile(dataDir,'anatomy','hughes');       % th
    fullfile(dataDir,'anatomy','vij');          % pv
    fullfile(dataDir,'anatomy','lim');          % sl
    fullfile(dataDir,'anatomy','veil');         % jv
    fullfile(dataDir,'anatomy','lopez');        % dl
    fullfile(dataDir,'anatomy','khazenzon');    % ak
    fullfile(dataDir,'anatomy','waskom');       % mw
    fullfile(dataDir,'anatomy','tiu');          % gt
    fullfile(dataDir,'anatomy','wsato');        % ws
    fullfile(dataDir,'anatomy','leung');        % ol
    fullfile(dataDir,'anatomy','lian');         % tl
    fullfile(dataDir,'anatomy','vitelli');      % mv
    fullfile(dataDir,'anatomy','martinez');     % vm
    fullfile(dataDir,'anatomy','bugno');        % ab
    fullfile(dataDir,'anatomy','wandell');      % bw
    fullfile(dataDir,'anatomy','brainard');     % dys_ab
    fullfile(dataDir,'anatomy','goodman');      % heb_ag
    fullfile(dataDir,'anatomy','Ayzenshtat');   % heb_aa
    fullfile(dataDir,'anatomy','Shambiro');     % heb_ls
    fullfile(dataDir,'anatomy','Toba');         % heb_toba
    fullfile(dataDir,'anatomy','Maya');         % heb_maya
    fullfile(dataDir,'anatomy','Yama');         % heb_yama
    fullfile(dataDir,'anatomy','Blta');         % heb_blta
    fullfile(dataDir,'anatomy','ben-shachar');  % heb_mibe
    fullfile(dataDir,'anatomy','Avbe');         % heb_avbe
    fullfile(dataDir,'anatomy','Nitr');         % heb_nitr
    fullfile(dataDir,'anatomy','Gilo');         % heb_gilo
    fullfile(dataDir,'anatomy','Eisa');         % heb_eisa
    fullfile(dataDir,'anatomy','Dael');         % heb_dael
    fullfile(dataDir,'anatomy','Taay');         % heb_taay
    fullfile(dataDir,'anatomy','Maaf');         % heb_maaf
    fullfile(dataDir,'anatomy','Brne');         % heb_brne
    fullfile(dataDir,'anatomy','Nagr');         % heb_nagr
    fullfile(dataDir,'anatomy','Avhi');         % heb_avhi
    fullfile(dataDir,'anatomy','Avar');         % heb_avar
    fullfile(dataDir,'anatomy','Nihe');         % heb_nihe
    fullfile(dataDir,'anatomy','Noco');         % heb_noco
    fullfile(dataDir,'anatomy','Shsh');         % heb_shsh
    };
% Copy anatomy once, we will sort it when we write it to FW
% src = oldanatomy;
% dst = fullfile(dataDir,'anatomy');
% copyfile(src,dst)

% The ROIs are in another location...Run the function with this oldDatDir
%{
ROIdir = '/share/wandell/data/anatomy';
for ii=2:length(bk.list_anatomy)
    dst  = fullfile(bk.list_anatomy{ii},'ROIs');
    fps  = split(dst,filesep);
    if ~exist(dst,'dir')
        mkdir(dst)
    end
    src = fullfile(ROIdir,fps{end-1},fps{end});
    if exist(src,'dir')
        disp(['Copying to ' dst ' ...'])
        srcs = dir(fullfile(src,'*'));
        for ns=3:length(srcs)
            s = fullfile(src,srcs(ns).name);
            copyfile(s,dst);
            disp(['... copied.\n'])
        end
    end
end
%}

% Now we need to edit the vANATOMYPATH and all the inplane directories etc
% in mrSESSION
if redefinePath;mrSESSION_redefine_paths(cr, bk);end




% directory with Hebrew, English, and Checkerboard retinotopy, all 7 degree
% in radius
bk.list_sessionHebrewRet_resize = {
    'do not have'   %     'jg'
    'do not have'   %     'ad'
    'do not have'   %     'cc'
    'do not have'   %     'jw'
    'do not have'   %     'rs'
    'do not have'   %     'sg'
    'do not have'   %     'th'
    'do not have'   %     'pv'
    'do not have'   %     'sl'
    'do not have'   %     'jv'
    'do not have'   %     'dl'
    'do not have'   %     'ak'
    'do not have'   %     'mw'
    'do not have'   %     'gt'
    'do not have'   %     'ws'
    'do not have'   %     'ol'
    'do not have'   %     'tl'
    'do not have'   %     'mv'
    'do not have'   %     'vm'
    'do not have'   %     'ab'
    'do not have'   %     'bw'
    'do not have'   %     'dys_ab'
    'do not have'   %     'heb_ag'
    'do not have'   %     'heb_aa'
    'do not have'   %     'heb_ls'
    'do not have'   %     'heb_toba'
    'do not have'   %     'heb_maya'
    'do not have'   %     'heb_yama'
    'do not have'   %     'heb_blta'
    'do not have'   %     'heb_mibe'
    'do not have'   %     'heb_avbe'
    fullfile(dataDir,'heb_pilot10','RetAndHebrewLoc_resize') %     'heb_nitr'
    'do not have'   %     'heb_gilo'
    'do not have'   %     'heb_eisa'
    fullfile(dataDir,'heb_pilot13','RetAndHebrewLoc_resize') %     'heb_dael'
    'do not have'   %     'heb_taay'
    fullfile(dataDir,'heb_pilot15','RetAndHebrewLoc_resize') %     'heb_maaf'
    'do not have'   %     'heb_brne'
    fullfile(dataDir,'heb_pilot17','RetAndHebrewLoc_resize') %     'heb_nagr'
    'do not have'   %     'heb_avhi'
    'do not have'   %     'heb_avar'
    'do not have'   %     'heb_nihe'
    'do not have'   %     'heb_noco'
    'do not have'   %     'heb_shsh'
    };
if checkCopy,checkAndCopy(bk.list_sessionHebrewRet_resize, bk.list_sessionHebrewRet_resize),end


% directory large and small faces and words. and in some cases also checkers
bk.list_sessionSizeRet = {
    fullfile(dataDir,'jg','20150701_wordEcc_retFaceWord');    % jg
    'do not have';                                                        % ad
    fullfile(dataDir,'cc','tiledLoc_sizeRet');                % cc
    fullfile(dataDir,'jw','20150701_wordEcc_retFaceWord');    % jw
    'do not have';                                               % rs
    'do not have';                                               % sg
    'do not have';                                               % th
    'do not have';                                               % pv
    'do not have';                                               % sl
    'do not have';                                               % jv
    'do not have';                                               % dl
    'do not have';                                               % ak
    fullfile(dataDir,'mw','tiledLoc_sizeRet');     % mw
    fullfile(dataDir,'gt','tiledLoc_sizeRet');     % gt
    fullfile(dataDir,'ws','tiledLoc_sizeRet');     % ws
    fullfile(dataDir,'ol','tiledLoc_sizeRet');     % ol
    fullfile(dataDir,'tl','Localizer_sizeRet');    % tl
    fullfile(dataDir,'mv','tiledLoc_sizeRet');     % mv
    fullfile(dataDir,'vm','tiledLoc_sizeRet');     % vm
    fullfile(dataDir,'ab','tiledLoc_sizeRet');     % ab
    fullfile(dataDir,'bw','tiledLoc_sizeRet');     % bw
    fullfile(dataDir,'dys_ab','tiledLoc_sizeRet'); % dys_ab
    'do not have';                                 % heb_ag
    'do not have';                                 % heb_aa
    'do not have';                                 % heb_ls
    'do not have';                                 % heb_toba
    'do not have';                                 % heb_maya
    'do not have';                                 % heb_yama
    'do not have';                                 % heb_blta
    'do not have';                                 % heb_mibe
    'do not have';                                 % heb_avbe
    'do not have';                                 % heb_nitr
    'do not have';                                 % heb_gilo
    'do not have';                                 % heb_eisa
    'do not have';                                 % heb_dael
    'do not have';                                 % heb_taay
    'do not have';                                 % heb_maaf
    'do not have';                                 % heb_brne
    'do not have';                                 % heb_nagr
    'do not have';                                 % heb_avhi
    'do not have';                                 % heb_avar
    'do not have';                                 % heb_nihe
    'do not have';                                 % heb_noco
    'do not have';                                 % heb_shsh
    };
if checkCopy,checkAndCopy(bk.list_sessionSizeRet, bk.list_sessionSizeRet),end


% directory with tiled localizer.
% directory large and small faces and words. and in some cases also checkers
bk.list_sessionTiledLoc = {
    fullfile(dataDir,'jg','tiledLoc');                        % jg
    'do not have';                                                          % ad
    fullfile(dataDir,'cc','tiledLoc_sizeRet');                % cc
    'do not have';                                                          % jw
    'do not have';                                                          % rs
    'do not have';                                                          % sg
    'do not have';                                                          % th
    'do not have';                                                          % pv
    'do not have';                                                          % sl
    'do not have';                                                          % jv
    'do not have';                                                          % dl
    'do not have'                                                           % ak
    fullfile(dataDir,'mw','Localizer_16Channel');             % mw
    fullfile(dataDir,'gt','Localizer_16Channel');             % gt
    fullfile(dataDir,'ws','tiledLoc_sizeRet');                % ws
    fullfile(dataDir,'ol','tiledLoc_sizeRet');                % ol
    fullfile(dataDir,'tl','Localizer_sizeRet');               % tl
    fullfile(dataDir,'mv','tiledLoc_sizeRet');                % mv
    fullfile(dataDir,'vm','tiledLoc_sizeRet');                % vm
    fullfile(dataDir,'ab','tiledLoc_sizeRet');                % ab
    fullfile(dataDir,'bw','tiledLoc_sizeRet');                % bw
    fullfile(dataDir,'dys_ab','tiledLoc_sizeRet');            % dys_ab
    fullfile(dataDir,'heb_pilot01','Analyze_pseudoInplane');  % heb_ag
    fullfile(dataDir,'heb_pilot02','RetAndLoc');              % heb_aa
    fullfile(dataDir,'heb_pilot03','RetAndLoc');              % heb_ls
    fullfile(dataDir,'heb_pilot04','RetAndLoc_noXform');      % heb_toba
    fullfile(dataDir,'heb_pilot05','RetAndLoc');              % heb_maya
    fullfile(dataDir,'heb_pilot06','writtenWordLoc_OTperscription');  % heb_yama
    fullfile(dataDir,'heb_pilot07','RetAndHebrewLoc');                % heb_blta
    'do not have';                                                  % heb_mibe
    'do not have';                                                  % heb_avbe
    'do not have';                                                  % heb_nitr
    'do not have';                                                  % heb_gilo
    'do not have';                                                  % heb_eisa
    'do not have';                                                  % heb_dael
    'do not have';                                                  % heb_taay
    'do not have';                                                  % heb_maaf
    'do not have';                                                  % heb_brne
    'do not have';                                                  % heb_nagr
    'do not have';                                                  % heb_avhi
    'do not have'   %     'heb_avar'
    'do not have'   %     'heb_nihe'
    'do not have'   %     'heb_noco'
    'do not have'   %     'heb_shsh'
    };
if checkCopy,checkAndCopy(bk.list_sessionTiledLoc, bk.list_sessionTiledLoc),end

% test retest reliability: checkers and words
bk.list_sessionTestRetest = {
    'do not have';     % jg
    'do not have';     % ad
    fullfile(dataDir,'cc','ret_testRetest');     % cc
    'do not have';     % jw
    'do not have';     % rs
    'do not have';     % sg
    'do not have';     % th
    'do not have';     % pv
    'do not have';     % sl
    'do not have';     % jv
    'do not have';     % dl
    'do not have';     % ak
    'do not have';     % mw
    'do not have';     % gt
    'do not have';     % ws
    'do not have';     % ol
    'do not have';     % tl
    'do not have';     % mv
    'do not have';     % vm
    fullfile(dataDir,'ab','ret_testRetest');     % ab
    'do not have';     % bw
    'do not have';     % dys_ab
    'do not have';     % heb_ag
    'do not have';     % heb_aa
    'do not have';     % heb_ls
    'do not have';     % heb_toba
    'do not have';     % heb_maya
    'do not have';     % heb_yama
    'do not have';     % heb_blta
    'do not have';     % heb_mibe
    'do not have';     % heb_avbe
    'do not have';     % heb_nitr
    'do not have';     % heb_gilo
    'do not have';     % heb_eisa
    'do not have';     % heb_dael
    'do not have';     % heb_taay
    'do not have';     % heb_maaf
    'do not have';     % heb_brne
    'do not have';     % heb_nagr
    'do not have';     % heb_avhi
    'do not have'   %     'heb_avar'
    'do not have'   %     'heb_nihe'
    'do not have'   %     'heb_noco'
    'do not have'   %     'heb_shsh'
    };
if checkCopy,checkAndCopy(bk.list_sessionTestRetest, bk.list_sessionTestRetest),end


% list freesurfer roi directory
bk.list_fsDir = {
    'still need to fill in'                                   % jg
    fullfile(dataDir,'anatomy','aaron3');         % ad
    fullfile(dataDir,'anatomy','camacho3');       % cc
    fullfile(dataDir,'anatomy','wexler');         % jw
    fullfile(dataDir,'anatomy','schneider');      % rs
    fullfile(dataDir,'anatomy','gleberman');      % sg
    fullfile(dataDir,'anatomy','hughes');         % th
    fullfile(dataDir,'anatomy','vij');            % pv
    fullfile(dataDir,'anatomy','lim');            % sl
    fullfile(dataDir,'anatomy','veil');           % jv
    fullfile(dataDir,'anatomy','lopez');          % dl
    fullfile(dataDir,'anatomy','khazenzon');      % ak
    'do not have';                                            % mw
    fullfile(dataDir,'anatomy','gtiu');           % gt
    fullfile(dataDir,'anatomy','wsato');          % ws
    fullfile(dataDir,'anatomy','leung');          % ol
    fullfile(dataDir,'anatomy','lian');           % tl
    fullfile(dataDir,'anatomy','vitelli');        % mv
    fullfile(dataDir,'anatomy','martinez');       % vm
    fullfile(dataDir,'anatomy','bugno');          % ab
    fullfile(dataDir,'anatomy','wandell');        % bw
    fullfile(dataDir,'anatomy','brainard');       % dys_ab
    fullfile(dataDir,'anatomy','goodman');        % heb_ag
    fullfile(dataDir,'anatomy','Ayzenshtat');     % heb_aa
    fullfile(dataDir,'anatomy','Shambiro');       % heb_ls
    fullfile(dataDir,'anatomy','Toba');           % heb_toba
    fullfile(dataDir,'anatomy','Maya');           % heb_maya
    fullfile(dataDir,'anatomy','Yama');           % heb_yama
    fullfile(dataDir,'anatomy','Blta');           % heb_blta
    fullfile(dataDir,'anatomy','ben-shachar');    % heb_mibe
    fullfile(dataDir,'anatomy','Avbe');           % heb_avbe
    fullfile(dataDir,'anatomy','Nitr');           % heb_nitr
    fullfile(dataDir,'anatomy','Gilo');           % heb_gilo
    fullfile(dataDir,'anatomy','Eisa');           % heb_eisa
    fullfile(dataDir,'anatomy','Dael');           % heb_dael
    fullfile(dataDir,'anatomy','Taay');           % heb_taay
    fullfile(dataDir,'anatomy','Maaf');           % heb_maaf
    fullfile(dataDir,'anatomy','Brne');           % heb_brne
    fullfile(dataDir,'anatomy','Nagr');           % heb_nagr
    fullfile(dataDir,'anatomy','Avhi');           % heb_avhi
    fullfile(dataDir,'anatomy','Avar');           % heb_avar
    fullfile(dataDir,'anatomy','Nihe');           % heb_nihe
    fullfile(dataDir,'anatomy','Noco');           % heb_noco
    fullfile(dataDir,'anatomy','Shsh');           % heb_shsh
    };
% if checkCopy,checkAndCopy(bk.list_fsDir, bk.list_fsDir),end


% directories with afq
% bk.list_sessionAfq = {
%     'not yet collected'                                           % jg
%     fullfile(dataDir,'ad','20150717_dti_qmri)')        % ad
%     'not yet collected'                                           % cc
%     fullfile(dataDir,'jw','20150426_dti_qmri')        % jw
%     fullfile(dataDir,'rs','20150510_dti_qmri')        % rs
%     fullfile(dataDir,'sg','20150517_dti_qmri')        % sg
%     fullfile(dataDir,'th','20150514_dti_qmri')        % th
%     fullfile(dataDir,'pv','20150510_dti_qmri')        % pv
%     fullfile(dataDir,'sl','20150509_dti_qmri')        % sl
%     fullfile(dataDir,'jv','20150512_dti_qmri')        % jv
%     'not yet collected'                                                 % dl
%     'not yet collected'                                                 % ak
%     'have diffusion data but afq not run yet'                           % mw
%     'have diffusion data but afq not run yet'                           % gt
%     'have diffusion data but afq not run yet'                           % ws
%     'have diffusion data but afq not run yet'                           % ol
%     'have diffusion data but afq not run yet'                           % tl
%     'have diffusion data but afq not run yet'                           % mv
%     'not yet collected'                                                 % vm
%     'not yet collected'                                                 % ab
%     'not yet collected'
%     'not yet collected'
%     'not yet collected'
%     'not yet collected'
%     };

%% retinotopy bookKeeping
eightROIs = {    'WangAtlas_V1v_left';
                 'WangAtlas_V2v_left';
                 'WangAtlas_V3v_left';
                 'WangAtlas_hV4_left';
                 'WangAtlas_VO1_left';
                 'lVOTRC';
                 'WangAtlas_IPS0';
                 'WangAtlas_IPS1'};
eightROIs = {};             
bk.list_ROIs = {
    eightROIs;      % jg
    eightROIs;      % ad
    eightROIs;      % cc
    eightROIs;      % jw
    eightROIs;      % rs
    eightROIs;      % sg
    eightROIs;      % th
    eightROIs;      % pv
    eightROIs;      % sl
    eightROIs;      % jv
    eightROIs;      % dl
    eightROIs;      % ak
    eightROIs;      % mw
    eightROIs;      % gt
    eightROIs;      % ws
    eightROIs;      % ol
    eightROIs;      % tl
    eightROIs;      % mv
    eightROIs;      % vm
    eightROIs;      % ab
    eightROIs;      % bw
    eightROIs;      % dys_ab
    eightROIs;      % heb_ag
    eightROIs;      % heb_aa
    eightROIs;      % heb_ls
    eightROIs;      % heb_toba
    eightROIs;      % heb_maya
    eightROIs;      % heb_yama
    eightROIs;      % heb_blta
    eightROIs;      % heb_mibe
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs;
    eightROIs
};

prfModel = {'css'};

bk.list_prfModels = {
    prfModel;      % jg
    prfModel;      % ad
    prfModel;      % cc
    prfModel;      % jw
    prfModel;      % rs
    prfModel;      % sg
    prfModel;      % th
    prfModel;      % pv
    prfModel;      % sl
    prfModel;      % jv
    prfModel;      % dl
    prfModel;      % ak
    prfModel;      % mw
    prfModel;      % gt
    prfModel;      % ws
    prfModel;      % ol
    prfModel;      % tl
    prfModel;      % mv
    prfModel;      % vm
    prfModel;      % ab
    prfModel;      % bw
    prfModel;      % dys_ab
    prfModel;      % heb_ag
    prfModel;      % heb_aa
    prfModel;      % heb_ls
    prfModel;      % heb_toba
    prfModel;      % heb_maya
    prfModel;      % heb_yama
    prfModel;      % heb_blta
    prfModel;      % heb_mibe
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel;
    prfModel
};



bk.list_rmName = [
    {{'Checkers';'Words';'FalseFont'}}      % jg
    {{'Checkers';'Words';'FalseFont'}}      % ad
    {{'Checkers';'Words';'FalseFont'}}      % cc
    {{'Checkers';'Words';'FalseFont'}}      % jw
    {{'Checkers';'Words';'FalseFont'}}      % rs
    {{'Checkers';'Words';'FalseFont'}}      % sg
    {{'Checkers';'Words';'FalseFont'}}      % th
    {{'Checkers';'Words';'FalseFont'}}      % pv
    {{'Checkers';'Words';'FalseFont'}}      % sl
    {{'Checkers';'Words';'FalseFont'}}      % jv
    {{'Checkers';'Words';'FalseFont'}}      % dl
    {{'Checkers';'Words';'FalseFont'}}      % ak
    {{'Checkers';'Words'}}      % mw 13
    {{'Checkers';'Words'}}      % gt
    {{'Checkers';'Words'}}      % ws
    {{'Checkers';'Words'}}      % ol
    {{'Checkers';'Words'}}      % tl 17
    {{'Checkers';'Words'}}      % mv
    {{'Checkers';'Words'}}      % vm
    {{'Checkers';'Words'}}      % ab
    {{'Checkers';'Words'}}      % bw
    {{'Checkers';'Words'}}      % dys_ab
    {{'Words_English';'Words_Hebrew'}}      % heb_ag
    {{'Words_English';'Words_Hebrew'}}      % heb_aa
    {{'Words_English';'Words_Hebrew'}}     % heb_ls
    {{'Words_English';'Words_Hebrew'}}      % heb_toba
    {{'Words_English';'Words_Hebrew'}}      % heb_maya
    {{'Words_English';'Words_Hebrew'}}      % heb_yama
    {{'Words_English';'Words_Hebrew'}}      % heb_blta
    {{'Words_English';'Words_Hebrew'}}      % heb_mibe
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    {{'Words_English';'Words_Hebrew'}}
    ];
%     {{'Checkers'}}      % heb_ag
%     {{'Checkers'}}      % heb_aa
%     {{'Checkers'}}     % heb_ls
%     {{'Checkers'}}      % heb_toba
%     {{'Checkers'}}      % heb_maya
%     {{'Checkers'}}      % heb_yama
%     {{'Checkers'}}      % heb_blta
%     {{'Checkers'}}      % heb_mibe
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     {{'Checkers'}}
%     ];
    
    
    
    
%     {{'Words_English';'Words_Hebrew'}}      % heb_ag
%     {{'Words_English';'Words_Hebrew'}}      % heb_aa
%     {{'Words_English';'Words_Hebrew'}}      % heb_ls
%     {{'Words_English';'Words_Hebrew'}}      % heb_toba
%     {{'Words_English';'Words_Hebrew'}}      % heb_maya
%     {{'Words_English';'Words_Hebrew'}}      % heb_yama
%     {{'Words_English';'Words_Hebrew'}}      % heb_blta
%     {{'Words_English';'Words_Hebrew'}}      % heb_mibe
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     {{'Words_English';'Words_Hebrew'}}
%     ];

bk.list_imFile_Checkers = {
    'Stimuli/images_8barswithblank_fliplr.mat';      % jg
    'Stimuli/images_8barswithblank_fliplr.mat';      % ad
    'Stimuli/images_8barswithblank_fliplr.mat';      % cc
    'Stimuli/images_8barswithblank_fliplr.mat';      % jw
    'Stimuli/images_8barswithblank_fliplr.mat';      % rs
    'Stimuli/images_8barswithblank_fliplr.mat';      % sg
    'Stimuli/images_8barswithblank_fliplr.mat';      % th
    'Stimuli/images_8barswithblank_fliplr.mat';      % pv
    'Stimuli/images_8barswithblank_fliplr.mat';      % sl
    'Stimuli/images_8barswithblank_fliplr.mat';      % jv
    'Stimuli/images_8barswithblank_fliplr.mat';      % dl
    'Stimuli/images_8barswithblank_fliplr.mat';      % ak
    'Stimuli/images_8barswithblank_fliplr.mat';      % mw
    'Stimuli/images_8barswithblank_fliplr.mat';      % gt
    'Stimuli/images_8barswithblank_fliplr.mat';      % ws
    'Stimuli/images_8barswithblank_fliplr.mat';      % ol
    'Stimuli/images_8barswithblank_fliplr.mat';      % tl
    'Stimuli/images_8barswithblank_fliplr.mat';      % mv
    'Stimuli/images_8barswithblank_fliplr.mat';      % vm
    'Stimuli/images_8barswithblank_fliplr.mat';      % ab
    'Stimuli/images_8barswithblank_fliplr.mat';      % bw
    'Stimuli/images_8barswithblank_fliplr.mat';      % dys_ab
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_ag
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_aa
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_ls
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_toba
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_maya
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_yama
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_blta
    'Stimuli/images_8barswithblank_fliplr.mat';      % heb_mibe
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    'Stimuli/images_8barswithblank_fliplr.mat';
    };

bk.list_imFile_KnK = {
    'Stimuli/images_knk_fliplr.mat';      % jg
    'Stimuli/images_knk_fliplr.mat';      % ad
    'Stimuli/images_knk_fliplr.mat';      % cc
    'Stimuli/images_knk_fliplr.mat';      % jw
    'Stimuli/images_knk_fliplr.mat';      % rs
    'Stimuli/images_knk_fliplr.mat';      % sg
    'Stimuli/images_knk_fliplr.mat';      % th
    'Stimuli/images_knk_fliplr.mat';      % pv
    'Stimuli/images_knk_fliplr.mat';      % sl
    'Stimuli/images_knk_fliplr.mat';      % jv
    'Stimuli/images_knk_fliplr.mat';      % dl
    'Stimuli/images_knk_fliplr.mat';      % ak
    'Stimuli/images_knk_fliplr.mat';      % mw
    'Stimuli/images_knk_fliplr.mat';      % gt
    'Stimuli/images_knk_fliplr.mat';      % ws
    'Stimuli/images_knk_fliplr.mat';      % ol
    'Stimuli/images_knk_fliplr.mat';      % tl
    'Stimuli/images_knk_fliplr.mat';      % mv
    'Stimuli/images_knk_fliplr.mat';      % vm
    'Stimuli/images_knk_fliplr.mat';      % ab
    'Stimuli/images_knk_fliplr.mat';      % bw
    'Stimuli/images_knk_fliplr.mat';      % dys_ab
    'Stimuli/images_knk_fliplr.mat';      % heb_ag
    'Stimuli/images_knk_fliplr.mat';      % heb_aa
    'Stimuli/images_knk_fliplr.mat';      % heb_ls
    'Stimuli/images_knk_fliplr.mat';      % heb_toba
    'Stimuli/images_knk_fliplr.mat';      % heb_maya
    'Stimuli/images_knk_fliplr.mat';      % heb_yama
    'Stimuli/images_knk_fliplr.mat';      % heb_blta
    'Stimuli/images_knk_fliplr.mat';      % heb_mibe
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    'Stimuli/images_knk_fliplr.mat';
    };


bk.list_paramsFile_KnK = {
    'Stimuli/params_knkfull_multibar_blank.mat';      % jg
    'Stimuli/params_knkfull_multibar_blank.mat';      % ad
    'Stimuli/params_knkfull_multibar_blank.mat';      % cc
    'Stimuli/params_knkfull_multibar_blank.mat';      % jw
    'Stimuli/params_knkfull_multibar_blank.mat';      % rs
    'Stimuli/params_knkfull_multibar_blank.mat';      % sg
    'Stimuli/params_knkfull_multibar_blank.mat';      % th
    'Stimuli/params_knkfull_multibar_blank.mat';      % pv
    'Stimuli/params_knkfull_multibar_blank.mat';      % sl
    'Stimuli/params_knkfull_multibar_blank.mat';      % jv
    'Stimuli/params_knkfull_multibar_blank.mat';      % dl
    'Stimuli/params_knkfull_multibar_blank.mat';      % ak
    'Stimuli/params_knkfull_multibar_blank.mat';      % mw
    'Stimuli/params_knkfull_multibar_blank.mat';      % gt
    'Stimuli/params_knkfull_multibar_blank.mat';      % ws
    'Stimuli/params_knkfull_multibar_blank.mat';      % ol
    'Stimuli/params_knkfull_multibar_blank.mat';      % tl
    'Stimuli/params_knkfull_multibar_blank.mat';      % mv
    'Stimuli/params_knkfull_multibar_blank.mat';      % vm
    'Stimuli/params_knkfull_multibar_blank.mat';      % ab
    'Stimuli/params_knkfull_multibar_blank.mat';      % bw
    'Stimuli/params_knkfull_multibar_blank.mat';      % dys_ab
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_ag
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_aa
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_ls
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_toba
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_maya
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_yama
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_blta
    'Stimuli/params_knkfull_multibar_blank.mat';      % heb_mibe
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    'Stimuli/params_knkfull_multibar_blank.mat';
    };

bk.list_paramsFile_Checkers = {
    'Stimuli/params_checkers.mat';      % jg
    'Stimuli/params_checkers.mat';      % ad
    'Stimuli/params_checkers.mat';      % cc
    'Stimuli/params_checkers.mat';      % jw
    'Stimuli/params_checkers.mat';      % rs
    'Stimuli/params_checkers.mat';      % sg
    'Stimuli/params_checkers.mat';      % th
    'Stimuli/params_checkers.mat';      % pv
    'Stimuli/params_checkers.mat';      % sl
    'Stimuli/params_checkers.mat';      % jv
    'Stimuli/params_checkers.mat';      % dl
    'Stimuli/params_checkers.mat';      % ak
    'Stimuli/params_checkers.mat';      % mw
    'Stimuli/params_checkers.mat';      % gt
    'Stimuli/params_checkers.mat';      % ws
    'Stimuli/params_checkers.mat';      % ol
    'Stimuli/params_checkers.mat';      % tl
    'Stimuli/params_checkers.mat';      % mv
    'Stimuli/params_checkers.mat';      % vm
    'Stimuli/params_checkers.mat';      % ab
    'Stimuli/params_checkers.mat';      % bw
    'Stimuli/params_checkers.mat';      % dys_ab
    'Stimuli/params_checkers.mat';      % heb_ag
    'Stimuli/params_checkers.mat';      % heb_aa
    'Stimuli/params_checkers.mat';      % heb_ls
    'Stimuli/params_checkers.mat';      % heb_toba
    'Stimuli/params_checkers.mat';      % heb_maya
    'Stimuli/params_checkers.mat';      % heb_yama
    'Stimuli/params_checkers.mat';      % heb_blta
    'Stimuli/params_checkers.mat';      % heb_mibe
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    'Stimuli/params_checkers.mat';
    };

% bk.list_paramsFile_Checkers = {
%     'Stimuli/20150113T200620.mat';      % jg
%     'Stimuli/20150113T200620.mat';      % ad
%     'Stimuli/20150113T200620.mat';      % cc
%     'Stimuli/20150113T200620.mat';      % jw
%     'Stimuli/20150113T200620.mat';      % rs
%     'Stimuli/20150113T200620.mat';      % sg
%     'Stimuli/20150113T200620.mat';      % th
%     'Stimuli/20150113T200620.mat';      % pv
%     'Stimuli/20150113T200620.mat';      % sl
%     'Stimuli/20150113T200620.mat';      % jv
%     'Stimuli/20150113T200620.mat';      % dl
%     'Stimuli/20150113T200620.mat';      % ak
%     'Stimuli/20150113T200620.mat';      % mw
%     'Stimuli/20150113T200620.mat';      % gt
%     'Stimuli/20150113T200620.mat';      % ws
%     'Stimuli/20150113T200620.mat';      % ol
%     'Stimuli/20150113T200620.mat';      % tl
%     'Stimuli/20150113T200620.mat';      % mv
%     'Stimuli/20150113T200620.mat';      % vm
%     'Stimuli/20150113T200620.mat';      % ab
%     'Stimuli/20150113T200620.mat';      % bw
%     'Stimuli/20150113T200620.mat';      % dys_ab
%     'Stimuli/20150113T200620.mat';      % heb_ag
%     'Stimuli/20150113T200620.mat';      % heb_aa
%     'Stimuli/20150113T200620.mat';      % heb_ls
%     'Stimuli/20150113T200620.mat';      % heb_toba
%     'Stimuli/20150113T200620.mat';      % heb_maya
%     'Stimuli/20150113T200620.mat';      % heb_yama
%     'Stimuli/20150113T200620.mat';      % heb_blta
%     'Stimuli/20150113T200620.mat';      % heb_mibe
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     'Stimuli/20150113T200620.mat';
%     };

bk.list_wSearch = [
    '3';      % jg
    '3';      % ad
    '3';      % cc
    '3';      % jw
    '3';      % rs
    '3';      % sg
    '3';      % th
    '3';      % pv
    '3';      % sl
    '3';      % jv
    '3';      % dl
    '3';      % ak
    '3';      % mw
    '3';      % gt
    '3';      % ws
    '3';      % ol
    '3';      % tl
    '3';      % mv
    '3';      % vm
    '3';      % ab
    '3';      % bw
    '3';      % dys_ab
    '3';      % heb_ag
    '3';      % heb_aa
    '3';      % heb_ls
    '3';      % heb_toba
    '3';      % heb_maya
    '3';      % heb_yama
    '3';      % heb_blta
    '3';      % heb_mibe
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    '3';
    ];


bk.list_stimSize = [
    16;      % jg
    16;      % ad
    16;      % cc
    16;      % jw
    16;      % rs
    16;      % sg
    16;      % th
    16;      % pv
    16;      % sl
    16;      % jv
    16;      % dl
    16;      % ak
    16;      % mw
    16;      % gt
    16;      % ws
    16;      % ol
    16;      % tl
    16;      % mv
    16;      % vm
    16;      % ab
    16;      % bw
    16;      % dys_ab
    7;      % heb_ag
    7;      % heb_aa
    7;      % heb_ls
    7;      % heb_toba
    7;      % heb_maya
    7;      % heb_yama
    7;      % heb_blta
    7;      % heb_mibe
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    7;
    ];
% a scan number in vista where the Knk scan was run (300 sec).
% IMPORTANT: this corresponds to list_sessionPath
bk.list_scanNum_Knk = [
    2;      % jg
    2;      % ad
    1;      % cc
    1;      % jw
    2;      % rs
    1;      % sg
    6;      % th
    1;      % pv
    1;      % sl
    1;      % jv
    1;      % dl
    4;      % ak
    0;      % mw
    0;      % gt
    0;      % ws
    0;      % ol
    0;      % tl
    0;      % mv
    0;      % vm
    0;      % ab
    0;      % bw
    0;      % dys_ab
    0;      % heb_ag
    0;      % heb_aa
    0;      % heb_ls
    0;      % heb_toba
    0;      % heb_maya
    0;      % heb_yama
    0;      % heb_blta
    0;      % heb_mibe
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    ];

% a scan number in vista where the Checkers scan was run (204 sec)
% IMPORTANT: this corresponds to list_sessionPath
bk.list_scanNum_Checkers = [
    1;      % jg
    1;      % ad
    2;      % cc
    2;      % jw
    1;      % rs
    2;      % sg
    7;      % th
    2;      % pv
    2;      % sl
    2;      % jv
    7;      % dl
    1;      % ak
    0;      % mw
    0;      % gt
    0;      % ws
    0;      % ol
    0;      % tl
    0;      % mv
    0;      % vm
    0;      % ab
    0;      % bw
    1;      % dys_ab
    0;      % heb_ag
    0;      % heb_aa
    0;      % heb_ls
    0;      % heb_toba
    0;      % heb_maya
    0;      % heb_yama
    0;      % heb_blta
    0;      % heb_mibe
    0;      % heb_avbe
    0;      % heb_nitr
    0;      % heb_gilo
    0;      % heb_eisa
    0;      % heb_dael
    0;      % heb_taay
    0;      % heb_maaf
    0;      % heb_brne
    0;      % heb_nagr
    0;      % heb_avhi
    0;      % heb_avar
    0;      % heb_nihe
    0;      % heb_noco
    0;      % heb_shsh
    ];

% a scan number in vista where the Knk scan was run (300 sec).
% IMPORTANT: this corresponds to list_sessionRet
bk.list_scanNum_Knk_sessionRet = [
    2;      % jg
    2;      % ad
    1;      % cc
    1;      % jw
    2;      % rs
    1;      % sg
    6;      % th
    1;      % pv
    1;      % sl
    1;      % jv
    1;      % dl
    4;      % ak
    2;      % mw
    3;      % gt
    3;      % ws
    3;      % ol
    1;      % tl
    3;      % mv
    3;      % vm
    3;      % ab
    1;      % bw
    3;      % dys_ab
    6;      % heb_ag
    6;      % heb_aa
    6;      % heb_ls
    6;      % heb_toba
    6;      % heb_maya
    0;      % heb_yama
    5;      % heb_blta
    5;      % heb_mibe
    5;      % heb_avbe
    5;      % heb_nitr
    5;      % heb_gilo
    5;      % heb_eisa
    5;      % heb_dael
    5;      % heb_taay
    5;      % heb_maaf
    5;      % heb_brne
    5;      % heb_nagr
    5;      % heb_avhi
    5;      % heb_avar
    5;      % heb_nihe
    5;      % heb_noco
    5;      % heb_shsh
    ];

% a scan number in vista where the checkers scan was run (204 sec).
% IMPORTANT: this corresponds to list_sessionRet
bk.list_scanNum_Checkers_sessionRet = [
    1;      % jg
    1;      % ad
    2;      % cc
    2;      % jw
    1;      % rs
    2;      % sg
    7;      % th
    2;      % pv
    2;      % sl
    2;      % jv
    7;      % dl
    1;      % ak
    1;      % mw
    1;      % gt
    1;      % ws
    1;      % ol
    1;      % tl
    1;      % mv
    1;      % vm
    1;      % ab
    1;      % bw
    1;      % dys_ab
    4;      % heb_ag
    4;      % heb_aa
    4;      % heb_ls
    4;      % heb_toba
    4;      % heb_maya
    0;      % heb_yama
    3;      % heb_blta
    3;      % heb_mibe
    3;      % heb_avbe
    3;      % heb_nitr
    3;      % heb_gilo
    3;      % heb_eisa
    3;      % heb_dael
    3;      % heb_taay
    3;      % heb_maaf
    3;      % heb_brne
    3;      % heb_nagr
    3;      % heb_avhi
    3;      % heb_avar
    3;      % heb_nihe
    3;      % heb_noco
    3;      % heb_shsh
    ];

% a scan number where the checkers scan was run (300 sec).
% IMPORTANT: this corresponds to list_sessionTestRetest
bk.list_scanNum_Checkers_sessionTestRetest = [
    0;
    0;
    1;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    1;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    ];

% a scan number where the knk scan was run (300 sec).
% IMPORTANT: this corresponds to list_sessionTestRetest
bk.list_scanNum_Knk_sessionTestRetest = [
    0;
    0;
    3;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    3;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    0;
    ];


%% colors
% colors!
bk.list_colors = {
    [0 1 .9];
    [1 0 .2];
    [.2 .9 .1];
    };

    [0.3922    0.3922    0.5020];  % purple
    [1.0000    0.5647    0.3922];  % orange
    [0    0.6353    0.3255];       % teal

%% standardized mesh names
% left mesh
'lh_inflated400_smooth1.mat';
% right mesh
'rh_inflated400_smooth1.mat';


%% general info, in one spot

% kgs localizer data (for seeing where vwfa and stuff is)
% (old one) /biac2/kgs/projects/Longitudinal/FMRI/Localizer/data/
% (after SNI move) /share/wandell/biac2/kgs/projects/Longitudinal/FMRI/Localizer/data

% kgs anatomy: for t1s, class files, rois and such
% /share/wandell/biac2/kgs/3Danat

% ROIs ==================================================================
% see bookKeeping_ROIs.m

% names of parameter maps STANDARD =======================================
% 'WordVAll.mat' <-------------
% 'WordNumberVAll.mat'
% 'FaceVAll.mat'
% 'varExp_CheckersMinusWords.mat'
% 'varExp_WordsMinusFalseFont.mat'
% 'varExp_WordExceeds.mat' -- ((words - checkers) + (words - falsefont)) / 2

% names of parameter maps TILED ==========================================

% WordVFace_Scrambled.mat <------------
% WordVFace.mat

% FaceVWord_Scrambled.mat <-----------
% FaceVScrambled.mat

% Proportion Variance Explained.mat
% Residual Variance.mat


% names of the mesh views ================================================
% ventral_lh,  ventral_rh
% lateral_lh, lateral_rh
% medial_lh, medial_rh
%
% Reading review . Sub 20. ----------------------------------------------
% Left lateral surface showing reading circuitry regions
% readingReview_leftLateral

% Left ventral surface showing VWFA
% readingReview_leftVentral


% namespaces =============================================================
% pp        -- preprocess
% mshimg    -- nathan's MeshImages code
% pmaps     -- parameter maps
% rm        -- things having to do with the rm model (i.e. loading into the view)

% dataType names and retModel names ======================================
% Checkers, Checkers-run1, Checkers-run2, Checkers-run3
% Words, Words-run1, Words-run2
% FalseFont, FalseFont-run1, FalseFont-run2, FalseFont-run3

% Ret model names are just "retModel-" with the name of the datatype
% appeneded after it, and the roi, if applicable. For example:

% retModel-Checkers
% retModel-Words
% retModel-FalseFont
% retModel-Combined - the average of all 3
%
% retModel-Words-run1: word stimulus, run 1, whole brain
% retModel-Words-run1-lh_VWFA_rl:  word stimulus, run 1, roi

%% Wang ROI names corresponding to number in the atlas
bk.list_wangRois = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};


%% colors
tmp                    = load(fullfile(defaultsDir,'list_colorsPerSub'));
bk.list_colorsPerSub   = tmp.list_colorsPerSub;
tmp                    = load(fullfile(defaultsDir,'list_colorsPerWangRois'));
bk.list_colorsWangRois = tmp.list_colorsWangRois;

%% Mesh views
% Reading Review --------------------------------------------------------
% Wang Atlas Figure -- Sub20
% WangVentralView
% WangDorsalView


% chdir(curdir)
end





%% Functions
function checkAndCopy(listCopy, listCheck)
dataDir    = '/black/localhome/glerma/TESTDATA/PRF-StimDependence/DATA';
oldbaseDir = '/biac4/wandell/data/';
olddataDir = fullfile(oldbaseDir, 'reading_prf');
oldanatomy = '/biac4/wandell/data/reading_prf/anatomy/';


if isequal(listCopy, listCheck)
    % Create a copy of the data that I can upload to FW later on
    for ii=1:length(listCopy)
        dst  = listCopy{ii};
        if strcmp(dst,'do not have') || strcmp(dst,'')
            % disp('do not have');
        else
            fps  = split(dst,filesep);
            subj = fullfile(dataDir,fps{end-1});
            if ~exist(subj,'dir')
                mkdir(subj)
            end
            if ~exist(dst,'dir')
                src = fullfile(olddataDir,fps{end-1},fps{end});
                if exist(src,'dir')
                    disp(['Copying to ' dst ' ...'])
                    copyfile(src,dst);
                    disp(['... copied.\n'])
                end
            end
        end
    end
else
    % Check if it is a different folder and if it is, copy it
    % Create a copy of the data that I can upload to FW later on
    if length(listCopy) ~= length(listCheck)
        error('Check different amount of folders, they should be the same')
    end
    for ii=1:length(listCopy)
        dst  = listCopy{ii};
        if strcmp(dst,'do not have') || strcmp(dst,'')
            % disp('do not have');
        else
            % Check that it was not copied before
            if ~strcmp(listCheck{ii}, dst)
                fps  = split(dst,filesep);
                subj = fullfile(dataDir,fps{end-1});
                if ~exist(subj,'dir')
                    mkdir(subj)
                end
                if ~exist(dst,'dir')
                    src = fullfile(olddataDir,fps{end-1},fps{end});
                    if exist(src,'dir')
                        disp(['Copying to ' dst ' ...'])
                        copyfile(src,dst);
                        disp(['... copied.\n'])
                    end
                end
            end
        end
    end
end
end
