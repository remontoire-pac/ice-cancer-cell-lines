function [] = iceopts(wf,wflag,uflag)
% function [] = iceopts(wf,wflag,uflag)
%  Set working folder for ICE and express options.
%          wf: working folder
%       wflag: set warnings on (default) or off
%       uflag: deploy ICE in 'use' mode (default) or 'build' mode
%   Created: PAC 2020-04-05
%   Last updated: PAC 2020-4-08 to edit documentation.
%
    
    if (nargin<3||isequal(uflag,[]))
        uflag = true;
    end
    assert(isscalar(uflag)&&islogical(uflag),'Input uflag must be a scalar logical.');
    if (nargin<2||isequal(wflag,[]))
        wflag = true;
    end
    assert(isscalar(wflag)&&islogical(wflag),'Input wflag must be a scalar logical.');
    if (nargin<1||isequal(wf,[]))
        wf = pwd();
    end
    cd(wf);
    if (wflag)
        warning('on','all');
    else
        warning('off','all');
    end
    if (uflag)
        icenv(wf,true);
    else
        icenv(wf,false);
    end
    
end
