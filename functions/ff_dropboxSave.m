function ff_dropboxSave(varargi)

%% Grabs the current figure and saves it to dropbox
% If there are no inputs, grabs the figure title and uses this as the
% filename
% Otherwise: ff_dropboxSave('title', 'theTitleNameThatYouWant')

%% If no inputs, grabs the current figure title and saves it to dropbox
% filename is title name
% saves both the .png and .fig file

dropboxDir = '~/Dropbox/TRANSFERIMAGES/';

%% input parser
p = inputParser; 
addOptional(p, 'title', get(get(gca, 'title'), 'String')); 
addOptional(p, 'saveto', './')); 
parse(p, varargin{:});
titleName = p.Results.title; 
saveto    = p.Results.saveto; 

%%

if iscell(titleName)
    titleName = ff_cellstring2string(titleName); 
end

saveas(gcf, fullfile(saveto, [titleName '.png']), 'png')
saveas(gcf, fullfile(saveto, [titleName '.fig']), 'fig')

end

