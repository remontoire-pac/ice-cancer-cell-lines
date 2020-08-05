%% setup workspace
wf = 'D:\job\ice'; addpath(wf); % path to iceopts.m, icenv.m, chkfile.m, tree.csv
iceopts(wf,false,false);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% find range of cell-line namespace from cell-line annotations file
assert(exist('ice000CellLineAnnotation.mat','file')==2, ...
             'Cell-line annotation MAT file required to fix column number.');
load ice000CellLineAnnotation.mat db*;
mxc = size(dbCellLineAnno,2);
clear db*;

%% load Gygi proteomics data (global MS) and extract headers
DM1 = chkfile('total-proteome-_v2-normalized-protein-abundance.csv','Gygi proteomics',',','string',true);
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
dmProteinMassSpec = nan(size(tm1,1),mxc);
dmProteinMassSpec(:,double(extractAfter(cn1u,'ACH-'))) = tm1;
clear c* t*;

%% process protein identifiers to a metadata table, index and sort data
mtProteinMassSpec = table('Size',[numel(rnC) 2],'VariableTypes',{'string','string'}, ...
                         'VariableNames',{'gene_symbol','uniprot_id'});
assert(isequal(unique(count(rnC,'_')),[2;3;4]),'Expecting exactly 2, 3, or 4 underscores.');
rnC = extractBefore(rnC,strlength(rnC));
assert(isequal(unique(count(rnC,'_')),[1;2;3]),'Expecting exactly 1, 2, or 3 underscores after extractBefore.');
rtC = count(rnC,'_'); rt1 = extractBefore(extractAfter(rnC,strlength(rnC)-7),2)=="_";
rnC(rtC==1) = regexprep(rnC(rtC==1),'_',"/",1); % replace only _ with / if 1 total
rnC(rtC==2&rt1) = regexprep(rnC(rtC==2&rt1),'_',"/",2); % replace second _ at -7 with / if 2 total 
rnC(rtC==2&~rt1) = regexprep(rnC(rtC==2&~rt1),'_',"/",1); % replace first _ not at -7 with / if 2 total
rnC(rtC==3) = regexprep(rnC(rtC==3),'_',"/",2); % replace middle _ with / if 3 total
rnC = replace(rnC,"_","-"); % replace all remaining _ with -
assert(isequal(unique(count(rnC,'/')),1),'Expecting exactly 1 slash after regexprep.');
assert(isequal(unique(count(rnC,'_')),0),'Expecting exactly 0 underscores after regexprep.');
assert(isequal(unique(count(rnC,'-')),[0;1;2]),'Expecting exactly 0, 1, or 2 dashes after regexprep.');
mtProteinMassSpec(:,:) = table(extractBefore(rnC,'/'),extractAfter(rnC,'/'));
mtProteinMassSpec.gene_symbol(mtProteinMassSpec.gene_symbol=="x") = missing();
[mtProteinMassSpec,so] = sortrows(mtProteinMassSpec,{'gene_symbol','uniprot_id'});
dmProteinMassSpec = single(dmProteinMassSpec(so,:));
mtProteinMassSpec.row_idx = (1:numel(rnC))';
mtProteinMassSpec = mtProteinMassSpec(:,[end 1:(end-1)]);
clear rnC rt* so;

%% write data and metadata variables to MAT file and CSV files
save build\mat\ice005GygiGlobalProteomics.mat *ProteinMassSpec;
mat2csv(5,dmProteinMassSpec,mtProteinMassSpec,'global_ms_proteomics',wf);
