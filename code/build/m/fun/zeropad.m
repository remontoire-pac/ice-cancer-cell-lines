function [sn] = zeropad(i,n)
% function [sn] = zeropad(i,n)
%  Make a non-negative integer i into a zero-padded string of length n.   
%   Created: PAC before 2010
%   Modified: PAC 2016-03-23
%   Modified: PAC 2020-04-06 to allow 0 -> '000' and deny abuse of n.
%   Last updated: PAC 2020-04-30 to prevent recast of logical i or n.
%  

    if (nargin<2)
        n = 3;
    end
    assert(~islogical(n)&&isscalar(n)&&isindex(n-1),'Input n must be an integer (min: n=2).');
    L = numel(num2str(i));
    assert(~islogical(i)&&isindex(i+1),'Input i must be a non-negative integer.');
    assert(L<=n,'Input i must be at most n digits (default: n=3).');
    if (L==n)
        sn = num2str(i);
    else
        sn = [repmat('0',1,(n-L)) num2str(i)];
    end
        
end
