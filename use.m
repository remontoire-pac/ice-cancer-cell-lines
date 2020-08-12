%% setup workspace
wf = 'D:\job\ice'; addpath(wf); % path to iceopts.m, icenv.m, chkfile.m, tree.csv
iceopts(wf,false,true);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% load environment from MAT files
clearvars -except a* wf;
awfb = dir('build\mat\ice*');
awnfb = numel(awfb);
awtfb = nan(1,numel(awfb));
for ati=1:awnfb
    tic;
    load(awfb(ati).name);
    awtfb(ati) = toc;
end
fprintf('Total time to load from MAT %i seconds.\n',round(sum(awtfb)));
clear ati *fb;
