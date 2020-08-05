function f = isindex(M)
% function [f] = isindex(M)
%  Test whether a matrix M contains only positive integers.
%   Created: PAC before 2009-05-13
%   Modified: PAC 2009-05-13
%   Modified: PAC 2020-04-06 to edit documentation.
%   Last updated: PAC 2020-04-30 to force M to be non-empty and finite.
%  

    assert(numel(M)>0,'Input M must have at least one element.');
    % explicit integers
    if (isinteger(M)&&(min(M(:))>=1)&&isfinite(M))
        f = true;
    % de-facto integers
    elseif (isnumeric(M)&&isequal(M,fix(M))&&(min(M(:))>=1)&&isfinite(M))
        f = true;
    else
        f = false;
    end
    
end
