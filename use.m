%% setup workspace
if (~exist('wf','var')) 
    wf = regexp(matlab.desktop.editor.getActiveFilename,filesep,'split');
    wf = strjoin(wf(1:(numel(wf)-1)),filesep); % ICE root folder
end
addpath(wf); % path to iceopts.m, icenv.m, chkfile.m, tree.csv
iceopts(wf,true,false);         % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% load environment from MAT files
clearvars -except a* wf;
awfb = dir(['build' filesep 'mat' filesep 'ice*']);
if isempty(awfb)
    error('Please place input mat files data in ./build/mat/ and run script again.');
else
    awnfb = numel(awfb);
    awtfb = nan(1,numel(awfb));
    for ati=1:awnfb
        tic;
        load(awfb(ati).name);
        awtfb(ati) = toc;
    end
    fprintf('Total time to load from MAT %i seconds.\n',round(sum(awtfb)));
    clear ati *fb;
end
