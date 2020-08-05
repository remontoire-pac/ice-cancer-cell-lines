function T = chkfile(fn,msg,del,tt,rvn)
% function T = chkfile(fn,msg,del,tt,rvn)
%  Verify file fn exists of with delimiter del, text type tt, headers rvn (T/F)
%        and read to table T, otherwise return message containing msg.
%   Created: PAC 2020-04-08
%   Modified: PAC 2020-04-09
%   Updated: PAC 2020-04-25
%   Last updated: PAC 2020-07-04 not compatible with MATLAB 2020a.
%        https://www.mathworks.com/help/matlab/ref/readtable.html
%


    if (nargin<5||isequal(rvn,[]))
        rvn = true;
    end
    assert(isscalar(rvn)&&islogical(rvn), ...
            'Input rv must be a logical scalar.');
    if (nargin<4||isequal(tt,[]))
        tt = 'string';
    end
    assert(isequal(tt,'string')||isequal(tt,'char'), ...
            'Input tt must be ''string'' or ''char''.');
    if (nargin<3||isequal(del,[]))
        del = ',';
    end
    assert(isequal(del,',')||isequal(del,'\t'), ...
            'Input del must be comma '','' or tab ''\t''.');
    if (nargin<2||isequal(msg,[]))
        msg = 'Necessary';
    end
    assert(ischar(msg)&&(size(msg,1)==1),'Label msg must be a character row vector.');
    assert(size(msg,2)<=namelengthmax,'Label msg must be shorter than %i characters.',1+namelengthmax);
    try
        T = readtable(fn,'FileType','text','Delimiter',del, ...
                        'TextType',tt,'ReadVariableNames',rvn);
    catch
        fprintf('%s input file %s not found. Exiting.\n',msg,fn);
        T = NaN;
    end
    
end
