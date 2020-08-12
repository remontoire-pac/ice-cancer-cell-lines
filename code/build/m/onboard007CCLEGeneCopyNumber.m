%% setup workspace
if (~exist('wf','var')) 
    wf = regexp(matlab.desktop.editor.getActiveFilename,filesep,'split');
    wf = strjoin(wf(1:(numel(wf)-4)),filesep); % ICE root folder
end
iceopts(wf,false,false);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% find range of cell-line namespace from cell-line annotations file
assert(exist('ice000CellLineAnnotation.mat','file')==2, ...
             'Cell-line annotation MAT file required to fix column number.');
load ice000CellLineAnnotation.mat db*;
mxc = size(dbCellLineAnno,2);
clear db*;

%% load CCLE gene-level copy-number data and extract headers
DM1 = chkfile('CCLE_gene_cn.csv','CCLE copy-number',',','string',true);
rnC = fields(DM1); rnC = string(rnC(2:(end-3)));
cn1 = DM1{:,DM1.Properties.VariableNames{1,1}}; cn1 = arxcellsyn(cn1); % cell lines are rows
ctm = numel(find(contains(cn1,'No CELL_SAMPLE_ID found')));
assert(ctm==0,'Achilles name not produced for %i input names.',ctm);
tm1 = (double(strrep(DM1{:,2:end},'NA','NaN')))'; % transpose
clear DM1;

%% process cell-line identifiers and shape data accordingly
[cn1s,cs] = sort(cn1); tm1 = tm1(:,cs); % initial column sort
[cn1u,cfo,cbi] = unique(cn1s); % find duplicate Achilles names
cld = unique([find(diff(cbi)==0) find(diff(cbi)==0)+1]);
if (numel(cld)>0)
    cldd = find([diff(cld)>1 1]); % identify sets of duplicates
    cldd = [1 cldd(1:(end-1))+1;cldd]; % find duplicate ranges
    for ci=1:size(cldd,2)
        tcld = cld(cldd(1,ci):cldd(2,ci));
        tm1(:,tcld(1)) = mean(tm1(:,tcld),2); % mean (for now)
    end
end
tm1 = tm1(:,cfo); % keep unique labels only
dmGeneCopyNumber = nan(size(tm1,1),mxc);
dmGeneCopyNumber(:,double(extractAfter(cn1u,'ACH-'))) = tm1;
clear c* t*;

%% process gene identifiers to a metadata table, index and sort data
mtGeneCopyNumber = table('Size',[numel(rnC) 2],'VariableTypes',{'string','doublenan'}, ...
                          'VariableNames',{'gene_or_locus_name','entrez_gene_id'});
rtC = count(rnC,'_');
assert(isequal(unique(count(rnC,'_')),[2;3;4]),'Expecting exactly 2, 3, or 4 underscores.');
rnC(rtC==4) = regexprep(rnC(rtC==4),'_',"-",2); % replace first _ with - only if 4 total
rnC(rtC==4) = regexprep(rnC(rtC==4),'_',"-",1); % replace second _ with - only if 4 total
rnC(rtC==3) = regexprep(rnC(rtC==3),'_',"-",1); % replace first _ with - only if 3 total
assert(isequal(unique(count(rnC,'_')),2),'Expecting exactly 2 underscores after regexprep.');
mtGeneCopyNumber(:,:) = table(extractBefore(rnC,'_'),extractBefore(extractAfter(rnC,'_'),'_'));
[mtGeneCopyNumber,so] = sortrows(mtGeneCopyNumber,{'gene_or_locus_name','entrez_gene_id'});
dmGeneCopyNumber = single(dmGeneCopyNumber(so,:));
mtGeneCopyNumber.row_idx = (1:numel(rnC))';
mtGeneCopyNumber = mtGeneCopyNumber(:,[end 1:(end-1)]);
assert(nnz(isnan(mtGeneCopyNumber.entrez_gene_id))==0, ...
                 'Missing gene IDs for some genes.');
clear r*C so;

%% write data and metadata variables to MAT file and CSV files
save build\mat\ice007CCLEGeneCopyNumber.mat *GeneCopyNumber;
mat2csv(7,dmGeneCopyNumber,mtGeneCopyNumber,'gene_level_copy_number',wf);
