%% setup workspace
wf = 'D:\job\ice'; addpath(wf); % path to iceopts.m, icenv.m, chkfile.m, tree.csv
iceopts(wf,false,false);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% run entire build-set with timing
awfc = dir('code\build\m\onboard*');
awnfc = numel(awfc);
awtfc = nan(1,numel(awfc));
for ati=1:awnfc
    tic;
    run(awfc(ati).name);
    awtfc(ati) = toc;
end
fprintf('Total time to build from src %i seconds.\n',round(sum(awtfc)));
clear ati *fc;

%% clear and rebuild from MAT with timing
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
