function M = castcol(tm,ncn,mxc)
% function M = castcol(tm,ncn,mxc)
%  Cast a compact matrix into specific columns of a sparser one with the same rows
%     (used to cast multiple different matrices into a common space of columns)
%         tm: M x P input matrix to be cast into M x N (N>=P) output matrix
%        ncn: positions columns of tm will occupy in output matrix
%        mxc: scalar maximum number of columns for output matrix
%   Created: PAC 2020-08-10
%   Last updated: PAC 2020-08-10
%

    assert(isscalar(mxc)&isindex(mxc),'Input mxc must be a scalar index.');
    assert(isvector(ncn)&isindex(ncn),'Input ncn must be a vector of indices.');
    assert(max(ncn)<=mxc,'Indices in input ncn must not exceed input mxc.');
    assert(isequal(sort(ncn),unique(ncn)),'Input ncn must contain distinct indices.');
    assert(ismatrix(tm)&&isnumeric(tm),'Input tm must be a numeric matrix.');
    assert(size(tm,2)==numel(ncn),'Input tm must have columns for each entry in input ncn.');
    M = nan(size(tm,1),mxc);
    M(:,ncn) = tm;
    
end
