function [] = mat2csv(ice,data,meta,stub,wf)
% function [] = mat2csv(ice,data,meta,stub,wf)
%  Write a standard output of matrix triples and row-label table to csv files.
%        ice: serial number for work product
%       data: logical or numeric data matrix
%       meta: table of metadata for the rows of data
%       stub: human-readable label shared by output files
%         wf: working folder
%   Created: PAC 2020-04-05
%   Modified: PAC 2020-04-08
%   Modified: PAC 2020-04-25 to for documentation and testing.
%   Modified: PAC 2020-04-30 to refuse logical ice, complex data; reorder asserts.
%   Last updated: PAC 2020-07-06 to write uint64 indices and single-precision float values.
%

    if (nargin<5||isequal(wf,[]))
        wf = pwd();
    end
    assert(exist(wf,'dir')==7,'%s must exist to root the _build_ path.',wf);
    cd(wf);
    icenv(wf,false);
    if (nargin<4||isequal(stub,[]))
        stub = 'cell_line_feature';
    end
    assert(ischar(stub)&&(size(stub,1)==1),'Label stub must be a character row vector.');
    assert(size(stub,2)<=namelengthmax,'Label stub must be shorter than %i characters.',1+namelengthmax);
    assert(~iskeyword(stub),'Label stub may not be a protected MATLAB keyword.');
    assert(isvarname(stub),'Label stub must start with a letter and contain only letters, digits, and underscores.');
    assert(istable(meta),'Input meta must be a table.');
    q = summary(meta);
    assert(numel(fields(q))>0,'Struct summary(meta) has no fields.');
    q = fields(q);
    q = q(1);
    assert(isequal(q,{'row_idx'}),'Input meta must start with a row_idx variable.');
    q = meta{:,1};
    assert(isnumeric(q)&&iscolumn(q),'The row_idx variable of meta must contain a numeric column vector.'); 
    assert(isequal(q,(1:numel(q))'),'The row_idx variable of meta must sequentially label its rows.');
    assert(ismatrix(data),'Input data must be a matrix.');
    assert(isreal(data)&&(islogical(data)||isnumeric(data)),'Input data must be logical or numeric.');
    [nr,nc] = size(data);
    assert(nr==size(meta,1),'Input matrix data and table meta must have matching row numbers.');
    assert(~islogical(ice)&&isindex(ice+1)&&isscalar(ice),'Input ice must be a scalar, non-negative integer.');
    ices = zeropad(ice,3);
    if (islogical(data))
        nTS = nnz(find(data));
        TS = table('Size',[nTS 2],'VariableTypes',{'uint64','uint64'}, ...
                                  'VariableNames',{'row_idx','col_idx'});
        [TS.row_idx,TS.col_idx] = find(data);
    else
        nTS = nnz(~isnan(data));
        TS = table('Size',[nTS 3],'VariableTypes',{'uint64','uint64','singlenan'}, ...
                                  'VariableNames',{'row_idx','col_idx','value'});
        [TS.row_idx,TS.col_idx] = find(~isnan(data));
        TS.value = data(~isnan(data));
    end
    writetable(meta,['build\csv\label\ice' ices '.rows.' stub '.' ...
               num2str(size(meta,1)) '.csv'],'FileType','text', ...
               'Delimiter',',','QuoteStrings',true,'WriteVariableNames',true);
    writetable(TS,['build\csv\triple\ice' ices '.' stub '_X_ccl.' ...
                   num2str(nr) '_X_' num2str(nc) '.csv'], ...
                   'FileType','text','Delimiter',',', ...
                   'QuoteStrings',false,'WriteVariableNames',true);

end
