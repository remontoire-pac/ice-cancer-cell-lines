%% setup workspace
if (~exist('wf','var')), wf = pwd(); end
assert(exist(wf,'dir')==7,'Directory %s not found.',wf); addpath(wf); 
iceopts(wf,false,false);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% prepare cell-line namespace using Project Achilles names (ACH-999999)
T = chkfile('cell-line-synon-names.csv','Cell-line synonym',',','string',true);
[tmcn,tmxc] = arxcellsyn(unique(T.CELL_SAMPLE_NAME(contains(T.CELL_NAME_TYPE,'Arx'))));

%% cross-check DepMap public file to ensure all names are present in namespace
Td = chkfile('sample_info.csv','CCLE/DepMap sample info',',','string',true);
assert(nnz(~ismember(arxcellsyn(unique(Td.DepMap_ID)),tmcn))==0,'Problem with DepMap indexing.');
assert(isequal(Td.DepMap_ID,unique(Td.DepMap_ID)),'Input DepMap indices not unique or not sorted.');
clear Td;

%% prepare tables of cell-line names and annotation metadata
mtCellLineAnno = chkfile('cell-line-anno-names.csv','CDDB annotation name',',','string',true);
assert(isequal((1:size(mtCellLineAnno,1))',mtCellLineAnno.annotation_id), ...
                'Incoming index does not match row_idx.');
mtCellLineAnno = addvars(mtCellLineAnno,mtCellLineAnno.annotation_id, ...
                         'Before',1,'NewVariableNames',{'row_idx'});
mtCellLineAnno = mtCellLineAnno(:,1:(end-1));
Tm = chkfile('map-cell-anno-names.csv','CDDB annotation mapping',',','string',true);
Tb = chkfile('cell-line-best-names.csv','CDDB preferred name',',','string',true);
assert(isequal(Tb.CELL_SAMPLE_ID,unique(Tb.CELL_SAMPLE_ID)),'Incoming CELL_SAMPLE_ID not unique.');

%% allocate and populate sparse matrix of annotation calls
dbCellLineAnno = logical(spalloc(size(mtCellLineAnno,1),tmxc,size(Tm,1)));
tct = double(extractAfter(tmcn,'ACH-'));
tmbn = tmcn;
tmid = nan(size(tmcn));
for ti=1:numel(tmcn)
    tmid(ti) = T.CELL_SAMPLE_ID(T.CELL_SAMPLE_NAME==tmcn(ti));
    tmbn(ti) = Tb.CELL_SAMPLE_NAME(Tb.CELL_SAMPLE_ID==tmid(ti));
    ty = tct(ti);
    tx = Tm.annotation_id(Tm.CELL_SAMPLE_ID==tmid(ti));
    dbCellLineAnno(tx,ty) = true;
end
clear ti tx ty T Tm Tb;

%% prepare mapping of Project Achilles names to highest-priority synonym
mtCellLineName = table('Size',[tmxc 4],'VariableTypes',{'uint64','doublenan','string','string'}, ...
                                'VariableNames',{'col_idx','master_ccl_id','achilles_name','priority_name'});
mtCellLineName.col_idx = (1:tmxc)';
mtCellLineName.master_ccl_id(tct) = tmid;
mtCellLineName.achilles_name(tct) = tmcn;
mtCellLineName.priority_name(tct) = tmbn;
clear tmxc tmid tmcn tmbn tct;

%% write data and metadata variables to MAT file and CSV files
save build\mat\ice000CellLineAnnotation.mat *CellLineAnno;
save build\mat\ice000CellLineAnnotation.mat -append mtCellLineName;
writetable(mtCellLineName,['build\csv\label\ice000.cols.cell_line_identity.' ...
           num2str(size(mtCellLineName,1)) '.csv'],'FileType','text', ...
           'Delimiter',',','QuoteStrings',true,'WriteVariableNames',true);
mat2csv(0,dbCellLineAnno,mtCellLineAnno,'cell_line_annotation',wf);
