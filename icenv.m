function [] = icenv(wf,uflag)
% function [] = icenv(wf,flag)
%  Setup working folders for ICE depending on 'use' or 'build' mode.
%          wf: working folder
%       uflag: deploy ICE in 'use' mode (default) or 'build' mode
%   Created: PAC 2020-04-04
%   Modified: PAC 2020-04-08
%   Last updated: PAC 2020-04-25 to use chkfile.m.
%

    if (nargin<2||isequal(uflag,[]))
        uflag = true;
    end
    assert(isscalar(uflag)&&islogical(uflag),'Input flag must be a scalar logical.');
    if (nargin<1||isequal(wf,[]))
        wf = pwd();
    end
    tree = chkfile('tree.csv','Path setup',',','char',true);
    assert(istable(tree),'chkfile.m could not find tree.csv.');
    cd(wf);
    if (uflag) 
        tree = tree(abs(tree.use) > 0,:);
        for ti=1:size(tree,1)
            tf = tree.folder{ti};
            tf = replace(tf,'\',filesep());
            tf = [wf filesep() tf];
            tn = tree.use(ti);
            if (tn<0)
                assert(exist(tf,'dir')==7,'%s required in _use_ mode.',tf);
            else
                if (exist(tf,'dir')==0)
                    fprintf('Creating %s for use mode.\n',tf);
                    mkdir(tf);
                end
            end
        end
    else
        tree = tree(abs(tree.build) > 0,:);
        for ti=1:size(tree,1)
            tf = tree.folder{ti};
            tf = replace(tf,'\',filesep());
            tf = [wf filesep() tf];
            tn = tree.build(ti);
            if (tn<0)
                assert(exist(tf,'dir')==7,'%s required in _build_ mode.',tf);
            else
                if (exist(tf,'dir')==0)
                    fprintf('Creating %s for build mode.\n',tf);
                    mkdir(tf);
                end
            end
        end
    end
    addpath(genpath(wf));
    
end
