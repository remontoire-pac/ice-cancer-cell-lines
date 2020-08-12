function M = spavg(dvec,ridx,cidx)
% function M = spavg(dvec,ridx,cidx)
%  Create a sparse matrix by averaging duplicate values indexed by vectors
%     (used to create a matrix from a table without assuming non-redundancy)
%       dvec: value vector with potentially redundant indices
%       ridx: index vector for output rows
%       cidx: index vector for output cols
%   Created: PAC 2020-08-10
%   Last updated: PAC 2020-08-10
%

    assert(iscolumn(dvec)&isnumeric(dvec),'Input dvec must be a numeric column vector.');
    assert(iscolumn(ridx)&isindex(ridx),'Input ridx must be a column vector of indices.');
    assert(iscolumn(cidx)&isindex(cidx),'Input cidx must be a column vector of indices.');
    assert(isequal(size(ridx),size(dvec)),'Inputs ridx and dvec must be the same size.');
    assert(isequal(size(cidx),size(dvec)),'Inputs cidx and dvec must be the same size.');
    tv1s = sparse(ridx,cidx,dvec,max(ridx),max(cidx));
    tv1n = sparse(ridx,cidx,1,max(ridx),max(cidx));
    M = full(tv1s./tv1n);
    
end
