function Tr = tabridx(T)
% function Tr = tabridx(T)
%  Add a row index column row_idx at the beginning of a table
%         T: input table without a row_idx column
%   Created: PAC 2020-08-10
%   Last updated: PAC 2020-08-10
%

    assert(istable(T),'Input T must be a table.');
    assert(max(contains(fields(T),'row_idx'))==0,'Input T must not contain a row_idx column.');
    T.row_idx = (1:size(T,1))';
    Tr = T(:,[end 1:(end-1)]);
    
end
